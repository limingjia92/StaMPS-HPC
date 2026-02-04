function ps_parms_default()
%PS_PARMS_DEFAULT (HPC Optimized Version)
%   Sets default parameters for StaMPS-HPC processing.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Framework Optimization: Minimized I/O overhead .
%   2. Parameter Cleanup: Removed obsolete legacy parameters
%   3. Annotation: Fully annotated every parameter for clarity and maintainability.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    parmfile_name = 'parms.mat';
    parmfile_path = parmfile_name;

    % --- 1. Load Parameters ---
    if exist(parmfile_name, 'file')
        parms = load(parmfile_name);
    elseif exist(['..' filesep parmfile_name], 'file')
        parmfile_path = ['..' filesep parmfile_name];
        parms = load(parmfile_path);
    else
        parms = struct('Created', date);
        parms.small_baseline_flag = 'n'; 
    end

    % Snapshot fields to track changes for logging
    parmfields_before = fieldnames(parms);
    num_fields = size(parmfields_before, 1);

    % =========================================================================
    % STEP 2: ESTIMATE PHASE NOISE
    % =========================================================================
    
    if ~isfield(parms,'restart_flag')
        parms.restart_flag = 0;             % Restart previous run for ps_est_gamma_quick
    end
    if ~isfield(parms,'max_topo_err')
        parms.max_topo_err = 20;            % Max uncorrelated DEM error (m). Pixels with error > this are dropped.
    end
    if ~isfield(parms,'filter_grid_size')
        parms.filter_grid_size = 50;        % Pixel size of grid (m) for resampling before filtering.
    end
    if ~isfield(parms,'filter_weighting')
        parms.filter_weighting = 'P-square';% Weighting scheme for filtering ('P-square' or 'SNR').
    end
    if ~isfield(parms,'clap_win')
        parms.clap_win = 32;                % CLAP filter window size. Together with grid size, determines area for phase est.
    end
    if ~isfield(parms,'clap_low_pass_wavelength')
        parms.clap_low_pass_wavelength = 800; % CLAP filter low-pass cut-off spatial wavelength (m).
    end
    if ~isfield(parms,'clap_alpha')
        parms.clap_alpha = 1;               % CLAP alpha: relative contribution of low-pass vs adaptive phase.
    end
    if ~isfield(parms,'clap_beta')
        parms.clap_beta = 0.3;              % CLAP beta term.
    end
    if ~isfield(parms,'gamma_change_convergence')
        parms.gamma_change_convergence = 0.005; % Convergence threshold for gamma estimation loop.
    end
    if ~isfield(parms,'gamma_max_iterations')
        parms.gamma_max_iterations = 3;     % Max iterations for gamma estimation.
    end

    % =========================================================================
    % STEP 3: PS SELECTION 
    % =========================================================================

    if ~isfield(parms,'reest_gamma_flag')
        % Recalculate phase/coh during selection. 'y'=accurate/slow, 'n'=inaccurate/fast.
        parms.reest_gamma_flag = 'y';
    end
    if ~isfield(parms,'use_fast_topofit')
        % Using Mex version topofit. 'y'=Mex/fast, 'n'=Matlab/slow.
        parms.use_fast_topofit = 'y';
    end
    if ~isfield(parms,'select_method')
        parms.select_method = 'DENSITY';    % Pixel selection method ('DENSITY' or 'PERCENT').
    end
    if ~isfield(parms,'density_rand')
        % Max acceptable spatial density (per km^2) of random-phase pixels.
        if strcmpi(parms.small_baseline_flag,'y')
            parms.density_rand = 2;  
        else
            parms.density_rand = 20; 
        end
    end
    if ~isfield(parms,'percent_rand')
        % Max acceptable percentage of selected pixels with random phase.
        if strcmpi(parms.small_baseline_flag,'y')
            parms.percent_rand = 1; 
        else
            parms.percent_rand = 20;
        end
    end
    if ~isfield(parms,'gamma_stdev_reject')
        parms.gamma_stdev_reject = 0;       % Std dev threshold for coherence estimation (0 = feature disabled).
    end
    if ~isfield(parms,'slc_osf')
        parms.slc_osf = 1;                  % Raw SLC oversampling factor.
    end

    % =========================================================================
    % STEP 4: PS WEEDING 
    % =========================================================================

    if ~isfield(parms,'weed_time_win')
        parms.weed_time_win = 730;          % Smoothing window (days) for estimating phase noise distribution.
    end
    if ~isfield(parms,'weed_max_noise')
        parms.weed_max_noise = inf;         % Max allowed noise threshold for a pixel (Max ifg noise).
    end
    if ~isfield(parms,'weed_standard_dev')
        % Threshold standard deviation. Pixels > this are dropped.
        if strcmpi(parms.small_baseline_flag,'y')
            parms.weed_standard_dev = 1.2;
        else 
            parms.weed_standard_dev = 1.0;
        end
    end
    if ~isfield(parms,'weed_zero_elevation')
        parms.weed_zero_elevation = 'n';    % Drop pixels with zero elevation (usually sea/invalid).
    end
    if ~isfield(parms,'weed_neighbours')
        parms.weed_neighbours = 'n';        % Spatially weed redundant neighbours, keeping highest coherence (CPU intensive).
    end
    if ~isfield(parms,'drop_ifg_index')
        parms.drop_ifg_index = [];          % List of interferogram indices to exclude/drop.
    end
    
    % =========================================================================
    % STEP 5: PHASE CORRECTION & MERGING 
    % =========================================================================

    if ~isfield(parms,'merge_resample_size')
        % Grid size (m) for resampling during merge. 0 = no resampling.
        if strcmpi(parms.small_baseline_flag,'y')
            parms.merge_resample_size = 100; 
        else
            parms.merge_resample_size = 0;   
        end
    end
    if ~isfield(parms,'merge_standard_dev')
        parms.merge_standard_dev = inf;     % Std dev threshold for dropping resampled pixels during merge.
    end

    % =========================================================================
    % STEP 6: PHASE UNWRAPPING 
    % =========================================================================

    if ~isfield(parms,'unwrap_method')
        % Unwrapping algorithm (e.g., '3D', '3D_QUICK').
        if strcmpi(parms.small_baseline_flag,'y')
            parms.unwrap_method = '3D_QUICK';  
        else
            parms.unwrap_method = '3D';
        end
    end
    if ~isfield(parms,'unwrap_prefilter_flag')
        parms.unwrap_prefilter_flag = 'y';  % Prefilter phase before unwrapping ('y' recommended).
    end
    if ~isfield(parms,'unwrap_grid_size')
        parms.unwrap_grid_size = 200;       % Resampling grid spacing (m) if prefilter is 'y'.
    end
    if ~isfield(parms,'unwrap_gold_n_win')
        parms.unwrap_gold_n_win = 32;       % Window size for Goldstein filter.
    end
    if ~isfield(parms,'unwrap_time_win')
        parms.unwrap_time_win = 730;        % Smoothing window (days) for estimating phase noise distribution.
    end
    if ~isfield(parms,'unwrap_gold_alpha')
        parms.unwrap_gold_alpha = 0.8;      % Alpha value for Goldstein filter.
    end
    if ~isfield(parms,'unwrap_patch_phase')
        parms.unwrap_patch_phase = 'n';     % Use patch phase from Step 3 as prefiltered phase ('n' recommended).
    end
    if ~isfield(parms,'unwrap_hold_good_values')
        parms.unwrap_hold_good_values = 'n';% SBAS: Hold good pixels from first unwrap using sb_identify_good_pixels.
    end
    if ~isfield(parms,'unwrap_la_error_flag')
        parms.unwrap_la_error_flag = 'y';   % Account for Look Angle error during unwrapping.
    end
    if ~isfield(parms,'unwrap_spatial_cost_func_flag')
        parms.unwrap_spatial_cost_func_flag = 'n'; % Enable spatial cost function constraint in 3D unwrapping.
    end
    if ~isfield(parms,'subtr_tropo')
        parms.subtr_tropo = 'n';            % Remove tropospheric estimate ('y' or 'n').
    end
    if ~isfield(parms,'tropo_method')
        parms.tropo_method = 'a_l';         % Method for troposphere estimate (e.g., 'a_l', 'a_e').
    end

    % =========================================================================
    % STEP 7: ESTIMATE SCLA ERROR
    % =========================================================================
    
    if ~isfield(parms,'scla_drop_index')
        parms.scla_drop_index = [];         % Interferograms listed here excluded from SCLA calculation.
    end
    if ~isfield(parms,'scla_deramp')
        parms.scla_deramp = 'n';            % Estimate phase ramp for each ifg ('y'/'n').
    end
    if ~isfield(parms,'scla_method')
        parms.scla_method = 'L2';           % Method for estimating SCLA (e.g., 'L2').
    end
    if ~isfield(parms,'scn_wavelength')
        parms.scn_wavelength = 100;         % Spatially correlated noise wavelength (m).
    end
    if ~isfield(parms,'scn_time_win')
        parms.scn_time_win = 365;           % Time window for SCN estimation (days).
    end
    if ~isfield(parms,'scn_deramp_ifg')
        parms.scn_deramp_ifg = [];          % Indices of specific interferograms to apply deramping.
    end
    if ~isfield(parms,'scn_kriging_flag')
        parms.scn_kriging_flag = 'n';       % Use kriging-based SCN filtering (computationally expensive).
    end

    % =========================================================================
    % PLOTTING & MISC
    % =========================================================================
    
    if ~isfield(parms,'ref_lon')
        parms.ref_lon = [-inf,inf];         % Reference area longitude range for mean velocity.
    end
    if ~isfield(parms,'ref_lat')
        parms.ref_lat = [-inf,inf];         % Reference area latitude range.
    end
    if ~isfield(parms,'ref_centre_lonlat')
        parms.ref_centre_lonlat = [0,0];    % Center coordinates [lon, lat] of reference area.
    end
    if ~isfield(parms,'ref_radius')
        parms.ref_radius = inf;             % Radius of reference area.
    end
    if ~isfield(parms,'plot_dem_posting')
        parms.plot_dem_posting = 90;        % Resampling resolution (m) for background DEM plotting.
    end
    if ~isfield(parms,'plot_scatterer_size')
        parms.plot_scatterer_size = 120;    % Scatterer size in map units (m).
    end
    if ~isfield(parms,'plot_pixels_scatterer')
        parms.plot_pixels_scatterer = 3;    % Screen pixels occupied by each scatterer symbol.
    end
    if ~isfield(parms,'plot_color_scheme')
        parms.plot_color_scheme = 'inflation'; % Colormap: 'inflation', 'deflation', 'gray', 'GMT_relief', etc.
    end
    if ~isfield(parms,'shade_rel_angle')
        parms.shade_rel_angle = [90,45];    % Light angle [azimuth, elevation] for shaded relief DEM.
    end
    if ~isfield(parms,'lonlat_offset')
        parms.lonlat_offset = [0,0];        % Manual [lon, lat] offset to align PS points with DEM.
    end

    % =========================================================================
    % EXTERNAL FILE HANDLING (HPC Optimized)
    % =========================================================================

    % 1. Lambda (Wavelength)
    lambdaname = 'lambda.1.in';
    if ~isfield(parms,'lambda')
        found_lambda = false;
        if exist(lambdaname, 'file')
            parms.lambda = load(lambdaname);
            found_lambda = true;
        elseif exist(['..' filesep lambdaname], 'file')
            parms.lambda = load(['..' filesep lambdaname]);
            found_lambda = true;
        end
        if ~found_lambda, parms.lambda = NaN; end
    end

    % 2. Heading
    headingname = 'heading.1.in';
    if ~isfield(parms,'heading')
        found_heading = false;
        if exist(headingname, 'file')
             parms.heading = load(headingname);
             found_heading = true;
        elseif exist(['..' filesep headingname], 'file')
             parms.heading = load(['..' filesep headingname]);
             found_heading = true;
        end
        if ~found_heading, parms.heading = NaN; end
    end

    % 3. Processor Type
    if ~isfield(parms, 'insar_processor')
        processor_file = 'processor.txt';
        proc_path = '';
        
        if exist(processor_file, 'file')
            proc_path = processor_file;
        elseif exist(['..' filesep processor_file], 'file')
            proc_path = ['..' filesep processor_file];
        end
        
        if isempty(proc_path)
            parms.insar_processor = 'isce'; % HPC Default
        else
            processor = fileread(proc_path);
            processor = strtrim(processor);
            valid_processors = {'gamma', 'isce', 'gmtsar'};
            if ~ismember(lower(processor), valid_processors)
                warning('Stamps:ProcessorType', 'Non-standard processor "%s" found. Proceeding.', processor);
            end
            parms.insar_processor = processor;
        end
    end

    % =========================================================================
    % SAVE CHANGES
    % =========================================================================
    
    parmfields = fieldnames(parms);
    % Check if any new fields were added
    if size(parmfields, 1) ~= num_fields
        try
            save(parmfile_path, '-struct', 'parms');
            
            % Optimized Logging
            log_buffer = {};
            for i = 1:size(parmfields, 1)
                if isempty(strmatch(parmfields{i}, parmfields_before)) 
                   parmname = parmfields{i};
                   val = parms.(parmname);
                   
                   val_str = '[]';
                   if isnumeric(val)
                       val_str = num2str(val);
                   elseif ischar(val)
                       val_str = val;
                   end
                   
                   log_buffer{end+1} = [parmname, ' = ', val_str];
                end
            end
            
            for k = 1:length(log_buffer)
                logit(log_buffer{k}, 0);
            end

        catch ME
            error('Stamps:WriteError', ...
                'CRITICAL: Could not update default parameters in %s.\nError: %s', parmfile_path, ME.message);
        end
    end
end