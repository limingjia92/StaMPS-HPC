function [ph_uw, msd] = uw_3d(tmp_ph_file, xy, day, ifgday_ix, bperp, options, unwrap_ifg_index)
%UW_3D Unwrap phase time series (HPC Block-IO Optimized Version).
%
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
% ======================================================================
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       2.0 (Block-IO Architecture)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. File-Based IPC: Replaced massive in-memory matrix passing with a 
%      disk-mapped pointer ('tmp_ph_file') to prevent parfor OOM crashes.
%   2. Vectorized Block-IO Mapping: Transformed the final Grid-to-PS mapping 
%      into a streamed Block-IO loop, strictly bounding peak memory usage.
%   3. Structural Consolidation: Unified single and multi-master workflows 
%      and merged 'uw_unwrap_from_grid' inline for cleaner execution.
%   4. Interface Standardization: Refactored inputs to use a unified 'options' 
%      struct, simplifying parameter passing and dependency management.
%
% ======================================================================
%   ORIGINAL HEADER & USAGE (StaMPS)
% ======================================================================
%   Original Author: Andy Hooper, June 2007
%
%   Usage:
%   [PH_UW, MSD] = UW_3D(TMP_PH_FILE, XY, DAY, IFGDAY_IX, BPERP, OPTIONS, UNWRAP_IFG_INDEX)
%
%   Inputs:
%   TMP_PH_FILE = String path to the mapped wrapped phase file (e.g., 'ph_w_tmp.mat'),
%                 replacing the legacy N x M in-memory matrix to save RAM.
%   XY  = N x 2 matrix of coordinates in metres
%        (optional extra column, in which case first column is ignored)
%   DAY = vector of image acquisition dates in days, relative to master
%   IFGDAY_IX = M x 2 matrix giving index to master and slave date in DAY
%        for each interferogram 
%   BPERP  = M x 1 vector giving perpendicular baselines 
%   OPTIONS = structure optionally containing following fields:
%      la_flag    = look angle error estimation flag (def='y')
%      scf_flag   = spatial cost function estimation flag (def='y')
%      master_day = serial date number of master (used for info msg only, def=0)
%      grid_size  = size of grid in m to resample data to (def=5)
%      prefilt_win = size of prefilter window in resampled grid cells (def=16)
%      time_win   = size of time filter window in days (def=365)
%      unwrap_method (def='3D' for single master, '3D_FULL' otherwise)
%      goldfilt_flag = Goldstein filtering, 'y' or 'n' (def='n')
%      gold_alpha    = alpha value for Goldstein filter (def=0.8)
%      lowfilt_flag  = Low pass filtering, 'y' or 'n' (def='n')
%      n_trial_wraps = no. phase cycles poss between neighbouring points due to la error (def=6)
%      temp          = optional M x 1 vector of temperature difference for each ifg (def=[])
%      n_temp_wraps  = no. phase cycles poss between neighbouring points due to temp diff (def=2)
%      variance      = N x 1 matrix of variance values           
%   UNWRAP_IFG_INDEX = 1D vector specifying the indices of IFGs to process.
%
%   Outputs:
%   PH_UW = unwrapped phase
%   MSD   = mean square deviation of phase
% ======================================================================

% =========================================================================
% 1. Initialization and Option Parsing
% =========================================================================
tic;

if nargin < 4
    error('Error: Not enough arguments. Usage: uw_3d(ph_file, xy, day, ifgday_ix, ...)')
end

if nargin < 5 || isempty(bperp), bperp = []; end
if nargin < 6, options = struct; end
if nargin < 7, unwrap_ifg_index = 1:size(ifgday_ix,1); end

% Set default options
defaults = struct(...
    'la_flag', 'y', 'scf_flag', 'y', 'master_day', 0, 'grid_size', 5, ...
    'prefilt_win', 16, 'time_win', 365, 'unwrap_method', '3D', ...
    'goldfilt_flag', 'n', 'gold_alpha', 0.8, 'lowfilt_flag', 'n', ...
    'n_trial_wraps', 6, 'temp', [], 'n_temp_wraps', 2, ...
    'max_bperp_for_temp_est', 100, 'variance', [], 'ph_uw_predef_file', []);

opt_fields = fieldnames(defaults);
for i = 1:length(opt_fields)
    field = opt_fields{i};
    if ~isfield(options, field), options.(field) = defaults.(field); end
end

if size(xy, 2) == 2, xy(:, 2:3) = xy(:, 1:2); end
if size(day, 1) == 1, day = day'; end
if strcmpi(options.unwrap_method, '3D'), options.lowfilt_flag = 'y'; end

% =========================================================================
% 2. Phase Resampling (Gridding)
% =========================================================================
% CRITICAL: Pass the file string and the valid indices to uw_grid_wrapped
uw_grid_wrapped(tmp_ph_file, xy, options, unwrap_ifg_index);

% =========================================================================
% 3. Topology Construction
% =========================================================================
uw_interp; 

% =========================================================================
% 4. Spatiotemporal Unwrapping (Core)
% =========================================================================
uw_sb_unwrap_space_time(day, ifgday_ix, bperp, options);

% =========================================================================
% 5. Statistical Cost Calculation
% =========================================================================
uw_stat_costs(options);

% =========================================================================
% 6. Final Unwrapping & Grid-to-PS Mapping (Block-IO Vectorized)
% =========================================================================
fprintf('Unwrapping from grid (Block-IO Streaming)...\n')

% Load grid index metadata
uw = load('uw_grid', 'nzix', 'n_ps', 'grid_ij');
uu = load('uw_phaseuw', 'ph_uw', 'msd');

n_ps_total = size(xy, 1);
n_ifg_proc = length(unwrap_ifg_index);

% Create lookup table: Grid Index -> PS Unwrapped Index
gridix = zeros(size(uw.nzix));
gridix(uw.nzix) = 1:uw.n_ps; 
grid_lin_idx = sub2ind(size(uw.nzix), uw.grid_ij(:,1), uw.grid_ij(:,2));
ix_vec_global = gridix(grid_lin_idx); 

% Allocate output variable (Safe size: ~30M * 240 IFGs = ~24GB in RAM)
ph_uw = nan(n_ps_total, n_ifg_proc, 'single');

% Setup Block-IO reading from the mapped phase file
m_tmp = matfile(tmp_ph_file);
row_block_size = 1000000;
n_blocks = ceil(n_ps_total / row_block_size);
has_predef = isfield(options, 'ph_uw_predef_file') && ~isempty(options.ph_uw_predef_file);

fprintf('   Mapping grid results back to %d PS points in %d blocks...\n', n_ps_total, n_blocks);

% Calculate absolute continuous bounding box for Matfile I/O
min_ifg_idx = min(unwrap_ifg_index);
max_ifg_idx = max(unwrap_ifg_index);
relative_ifg_index = unwrap_ifg_index - min_ifg_idx + 1;

for b = 1:n_blocks
    r_start = (b-1)*row_block_size + 1;
    r_end = min(b*row_block_size, n_ps_total);
    rows = r_start:r_end;
    
    ix_vec_chunk = ix_vec_global(rows);
    valid_mask = (ix_vec_chunk > 0);
    
    if any(valid_mask)
        % 1. Load exact continuous data block from disk
        ph_raw_chunk_all = m_tmp.ph_w(rows, min_ifg_idx:max_ifg_idx);
        % Slice the non-contiguous IFGs in RAM
        ph_raw_chunk = ph_raw_chunk_all(:, relative_ifg_index);
        clear ph_raw_chunk_all;

        ph_smooth_chunk = uu.ph_uw(ix_vec_chunk(valid_mask), :);
        ph_raw_valid = ph_raw_chunk(valid_mask, :);
        
        % 2. Calculate residual and combine
        ph_resid = ph_raw_valid .* exp(-1i * ph_smooth_chunk);
        ph_uw_chunk_part = ph_smooth_chunk + angle(ph_resid);
        
        % 3. Apply predefined Phase Corrections (Cycle Skipping)
        if has_predef
            % Same proxy logic for predef phase
            predef_chunk_all = m_tmp.ph_uw_predef(rows, min_ifg_idx:max_ifg_idx);
            predef_chunk = predef_chunk_all(:, relative_ifg_index);
            clear predef_chunk_all;

            predef_valid = predef_chunk(valid_mask, :);
            predef_ix = ~isnan(predef_valid);
            
            if any(predef_ix(:))
                diff_val = ph_uw_chunk_part - predef_valid;
                meandiff = mean(diff_val(:), 'omitnan');
                cycle_shift = 2 * pi * round(meandiff / (2*pi));
                
                predef_valid = predef_valid + cycle_shift;
                ph_uw_chunk_part(predef_ix) = predef_valid(predef_ix);
            end
        end
        
        % 4. Insert back into pre-allocated output chunk
        chunk_out = nan(length(rows), n_ifg_proc, 'single');
        chunk_out(valid_mask, :) = ph_uw_chunk_part;
        ph_uw(rows, :) = chunk_out;
    end
end

msd = uu.msd;
fprintf('HPC Processing Complete. Total time: %.2f s\n', toc);

end