function stamps(start_step, end_step, patch_list_file)
%STAMPS Main control function for StaMPS-HPC processing chain.
%   
%   Usage:
%   stamps(1, 8)                 -> Standard full run (Auto-detects mode)
%   stamps(1, 5)                 -> Run Steps 1-5 (Single patch or Master dir)
%   stamps(5.5, 5.5)             -> Run Merge Step ONLY (ps_merge_patches + calc_ifg_std).
%   stamps(6, 8)                 -> Run Steps 6-8 (Post-processing/Merge)
%   stamps(1, 5, 'split_1.list') -> HPC Mode: Process specific patches in list
%
%   Inputs:
%   start_step      [Integer] Starting processing step (1-8). Default: 1.
%   end_step        [Integer] Ending processing step (1-8). Default: 8.
%   A subset of steps may be selected with START_STEP and/or END_STEP
%   STEP 1 = Initial load of data
%   STEP 2 = Estimate gamma 
%   STEP 3 = Select PS pixels
%   STEP 4 = Weed out adjacent pixels
%   STEP 5 = Correct wrapped phase for spatially-uncorrelated look angle error and merge patches
%   STEP 5.5 = Merge patches and calculate ifg standard deviation
%   STEP 6 = Unwrap phase
%   STEP 7 = Calculate spatially correlated look angle (DEM) error 
%   STEP 8 = Filter spatially correlated noise 
%
%   patch_list_file [String] (Optional) Path to a custom patch list file.
%                   - If specified: The function will ONLY process the patches listed 
%                   - If empty: The function auto-detects the mode:
%                     1. Master Directory: Uses 'patch.list' if found.
%                     2. Single Patch: Processes the current directory if no list exists.
%
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Workflow Restructuring: Decouples distributed patch processing to support parallel job array
%   2. Argument Simplification: Removed interactive flags in favor of auto-detection.
%   3. Logic Optimization: Enforcing optimized estimation algorithms and integrating internal loaders.
%   4. Batch Control: Support for external patch list files to enable precise task.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    % =========================================================================
    % 1. Initialization & Parameter Setup
    % =========================================================================
    
    nfill = 40;
    fillstr = [repmat('#', 1, nfill), '\n'];
    logit(fillstr);
    logit(' StaMPS/MTI Optimized HPC Version ');
    logit(fillstr);

    % Default Arguments
    if nargin < 1 || isempty(start_step), start_step = 1; end
    if nargin < 2 || isempty(end_step), end_step = 8; end
    if nargin < 3, patch_list_file = ''; end 

    % Load Parameters
    small_baseline_flag = getparm('small_baseline_flag');
    insar_processor     = getparm('insar_processor');
    scn_kriging_flag    = getparm('scn_kriging_flag');
    reest_gamma_flag    = getparm('reest_gamma_flag'); % Corrected name

    currdir = pwd;

    % =========================================================================
    % 2. Intelligent Mode Detection
    % =========================================================================
    is_multi_patch = false;
    
    if ~isempty(patch_list_file)
        % HPC Mode: Explicit patch list provided
        if exist(patch_list_file, 'file')
            is_multi_patch = true;
        else
            error('Provided patch list file not found: %s', patch_list_file);
        end
    else
        % Auto-detect Mode: Check for patch structure in current dir
        if exist('patch.list', 'file') || ~isempty(dir('PATCH_*'))
            is_multi_patch = true;
            patch_list_file = 'patch.list'; 
        else
            % Single Patch Mode (inside a patch folder)
            is_multi_patch = false;
        end
    end

    if is_multi_patch
        logit(sprintf('Mode: Multi-Patch / Batch Processing (List: %s)', patch_list_file));
    else
        logit('Mode: Single-Patch / Local Processing');
    end

    % =========================================================================
    % 3. PART 1: Per-Patch Processing (Steps 1 - 5a)
    % =========================================================================
    if start_step <= 5 && floor(start_step) == start_step
        
        % Determine directories to process
        if is_multi_patch
            patch_dirs = parse_patch_list(patch_list_file);
        else
            patch_dirs = {'.'}; 
        end

        for i = 1:length(patch_dirs)
            target_dir = patch_dirs{i};
            
            if ~strcmp(target_dir, '.')
                cd(target_dir);
            end
            
            patchsplit = strsplit(pwd, filesep);
            current_patch_name = patchsplit{end};
            
            % Initialize processing state tracker
            if exist('no_ps_info.mat', 'file') ~= 2
                stamps_step_no_ps = zeros([5 1]); 
                save('no_ps_info.mat', 'stamps_step_no_ps');
            else
                load('no_ps_info.mat'); 
            end

            % >>> STEP 1: Load Data <<<
            if start_step <= 1
                logit_step(1, current_patch_name);
                
                % Optimized Loading: Removed external data_inc logic
                if strcmpi(small_baseline_flag, 'y')
                    if strcmpi(insar_processor, 'gamma') || strcmpi(insar_processor, 'snap')
                        sb_load_initial_gamma;
                    elseif strcmpi(insar_processor, 'isce')
                         sb_load_initial_isce; 
                    elseif strcmpi(insar_processor, 'gsar')
                         sb_load_initial_gsar;
                    else
                        sb_load_initial;
                    end
                else
                    if strcmpi(insar_processor, 'gamma') || strcmpi(insar_processor, 'snap')
                        ps_load_initial_gamma;
                    elseif strcmpi(insar_processor, 'isce')
                        ps_load_initial_isce; 
                    elseif strcmpi(insar_processor, 'gsar')
                        ps_load_initial_gsar;
                    else
                        ps_load_initial;
                    end
                end
                
                stamps_step_no_ps(1:end) = 0;
                save('no_ps_info.mat', 'stamps_step_no_ps');
            end

            % >>> STEP 2: Estimate Gamma <<<
            if start_step <= 2 && end_step >= 2
                logit_step(2, current_patch_name);
                load('no_ps_info.mat');
                
                if stamps_step_no_ps(1) == 0
                    ps_est_gamma_quick; % Always use quick estimation
                else
                    logit('Skipping Step 2 (No PS from Step 1)');
                    stamps_step_no_ps(2) = 1;
                end
                save('no_ps_info.mat', 'stamps_step_no_ps');
            end

            % >>> STEP 3: PS Selection <<<
            if start_step <= 3 && end_step >= 3
                logit_step(3, current_patch_name);
                load('no_ps_info.mat');
                
                if stamps_step_no_ps(2) == 0
                    if strcmpi(reest_gamma_flag, 'y')
                        ps_select;
                    else
                        ps_select(1);
                    end
                else
                    logit('Skipping Step 3 (No PS from Step 2)');
                    stamps_step_no_ps(3) = 1;
                end
                save('no_ps_info.mat', 'stamps_step_no_ps');
            end

            % >>> STEP 4: Weeding <<<
            if start_step <= 4 && end_step >= 4
                logit_step(4, current_patch_name);
                load('no_ps_info.mat');
                
                if stamps_step_no_ps(3) == 0
                    if strcmpi(small_baseline_flag, 'y')
                        ps_weed(0, 1);
                    else
                        ps_weed;
                    end
                else
                    logit('Skipping Step 4 (No PS from Step 3)');
                    stamps_step_no_ps(4) = 1;
                end
                save('no_ps_info.mat', 'stamps_step_no_ps');
            end

            % >>> STEP 5a: Phase Correction <<<
            if start_step <= 5 && end_step >= 5
                logit_step(5, current_patch_name);
                load('no_ps_info.mat');
                
                if stamps_step_no_ps(4) == 0
                    ps_correct_phase;
                else
                    logit('Skipping Step 5 (No PS from Step 4)');
                    stamps_step_no_ps(5) = 1;
                end
                save('no_ps_info.mat', 'stamps_step_no_ps');
            end
            
            cd(currdir);
        end
    end

    % =========================================================================
    % 4. PART 2: Merge & Post-Processing (Steps 5b - 8)
    % =========================================================================
    
    explicit_merge = (start_step == 5.5 && end_step == 5.5);
    continuous_run = (start_step <= 5 && end_step >= 6);
    
    % --- Step 5.5: Merge / Transition ---
    if explicit_merge || (continuous_run && is_multi_patch)
         logit_step(5.5, 'Merging Patches');
         ps_merge_patches(2, patch_list_file);
         ps_calc_ifg_std;
         
         if explicit_merge
             return; % Exit if only asking for merge
         end
         
    elseif continuous_run && ~is_multi_patch
         ps_calc_ifg_std;
    end

    % --- Step 6-8: Post-Processing ---
    if end_step >= 6        
        logit('Starting Part 2: Post-Processing...');

        % >>> STEP 6: Unwrap <<<
        if start_step <= 6 && end_step >= 6
            logit_step(6, 'Unwrapping Phase');
            ps_unwrap;
            if strcmpi(small_baseline_flag, 'y')
                sb_invert_uw;
            end
        end

        % >>> STEP 7: SCLA Error <<<
        if start_step <= 7 && end_step >= 7
            logit_step(7, 'SCLA Error Estimation');
            if strcmpi(small_baseline_flag, 'y')
                ps_calc_scla(1, 1);   
                ps_smooth_scla(1);
                ps_calc_scla(0, 1);   
            else
                ps_calc_scla(0, 1);
                ps_smooth_scla;
            end
        end

        % >>> STEP 8: Filtering <<<
        if start_step <= 8 && end_step >= 8
            logit_step(8, 'Filtering Noise');
            if strcmpi(scn_kriging_flag, 'y')
                ps_scn_filt_krig;
            else
                ps_scn_filt;
            end
        end
    end
    
    logit('StaMPS Processing Complete.');
end

% =========================================================================
% Helper Functions
% =========================================================================

function logit_step(step_num, Text)
    fprintf('\n########################################\n');
    fprintf(' StaMPS Step %g: Processing %s \n', step_num, Text);
    fprintf('########################################\n');
end

function list = parse_patch_list(filename)
    % Parse patch list file into cell array
    fid = fopen(filename);
    if fid == -1
        error('Cannot open patch list file: %s', filename);
    end
    data = textscan(fid, '%s');
    list = data{1};
    fclose(fid);
end