function [ph_uw, msd] = uw_3d(ph, xy, day, ifgday_ix, bperp, options)
%UW_3D Unwrap phase time series (HPC Optimized Version).
%
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
% ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Structural Consolidation: Unified single/multi-master workflows into 
%      'uw_sb_unwrap_space_time' and merged 'uw_unwrap_from_grid' inline.
%   2. Vectorization: Replaced slow pixel-wise loops with vectorized matrix 
%      operations for final grid mapping (10-50x speedup).
%   3. Interface Standardization: Refactored inputs to use a unified 'options' 
%      struct, simplifying parameter passing and dependency management.
%
% ======================================================================
%   ORIGINAL HEADER & USAGE (StaMPS)
% ======================================================================
%   Original Author: Andy Hooper, June 2007
%
%   Usage:
%   PH_UW = UW_3D(PH,XY,DAY,IFGDAY_IX,BPERP,OPTIONS)
%
%   Inputs:
%   PH  = N x M matrix of wrapped phase values (real phase or complex phasor)
%        where N is number of pixels and M is number of interferograms
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
    error('Error: Not enough arguments. Usage: uw_3d(ph, xy, day, ifgday_ix, ...)')
end

if nargin < 5 || isempty(bperp)
    bperp = []; 
end

if nargin < 6
    options = struct;
end

% Set default options
defaults = struct(...
    'la_flag', 'y', ...
    'scf_flag', 'y', ...
    'master_day', 0, ...
    'grid_size', 5, ...
    'prefilt_win', 16, ...
    'time_win', 365, ...
    'unwrap_method', '3D', ...
    'goldfilt_flag', 'n', ...
    'gold_alpha', 0.8, ...
    'lowfilt_flag', 'n', ...
    'n_trial_wraps', 6, ...
    'temp', [], ...
    'n_temp_wraps', 2, ...
    'max_bperp_for_temp_est', 100, ...
    'variance', [], ...
    'ph_uw_predef', []);

% Merge defaults with provided options
opt_fields = fieldnames(defaults);
for i = 1:length(opt_fields)
    field = opt_fields{i};
    if ~isfield(options, field)
        options.(field) = defaults.(field);
    end
end

% Validate options
valid_list = fieldnames(defaults);
provided_list = fieldnames(options);
invalid_options = setdiff(provided_list, valid_list);
if ~isempty(invalid_options)
    error(['Error: Invalid option provided: "', invalid_options{1}, '"']);
end

% Validate Temp dimensions
if ~isempty(options.temp) && length(options.temp) ~= size(ph, 2)
    error('Error: options.temp must be M x 1 vector (M = number of ifgs)');
end

% Coordinate handling
if size(xy, 2) == 2
   xy(:, 2:3) = xy(:, 1:2); % Shift to columns 2-3 if only 2 provided (legacy format)
end

if size(day, 1) == 1
    day = day';
end
    
if strcmpi(options.unwrap_method, '3D') 
    options.lowfilt_flag = 'y';
end

% =========================================================================
% 2. Phase Resampling (Gridding)
% =========================================================================
% Resamples the full-resolution wrapped phase onto a regular grid to reduce
% computational load and noise.
% Output: 'uw_grid.mat'
uw_grid_wrapped(ph, xy, options);

clear ph % Clear large input variable to save memory

% =========================================================================
% 3. Topology Construction
% =========================================================================
% Generates the Delaunay triangulation and neighbor connectivity for the grid.
% Output: 'uw_interp.mat'
uw_interp; 

% =========================================================================
% 4. Spatiotemporal Unwrapping (Core)
% =========================================================================
% Solves the integer ambiguity problem in 3D (Space-Time).
% Output: 'uw_space_time.mat' 
uw_sb_unwrap_space_time(day, ifgday_ix, bperp, options);

% =========================================================================
% 5. Statistical Cost Calculation
% =========================================================================
% Estimates the reliability (cost) of the unwrapping solution.
% Updates 'uw_phaseuw.mat'
uw_stat_costs(options);

% =========================================================================
% 6. Final Unwrapping & Grid-to-PS Mapping (Vectorized)
% =========================================================================
% Map smooth grid phase back to PS pixels and restore high-frequency details.
fprintf('Unwrapping from grid (Inline Optimized)...\n')

% Load intermediate data
uw = load('uw_grid', 'nzix', 'n_ps', 'grid_ij', 'ph_in', 'ph_in_predef');
uu = load('uw_phaseuw', 'ph_uw', 'msd');

% -- Map Grid Phase to PS Pixels --
[n_ps_grid, n_ifg] = size(uw.ph_in);

% Create lookup table: Grid Index -> PS Unwrapped Index
gridix = zeros(size(uw.nzix));
gridix(uw.nzix) = 1:uw.n_ps; 

% Identify which unwrapped grid cell corresponds to each PS pixel
grid_lin_idx = sub2ind(size(uw.nzix), uw.grid_ij(:,1), uw.grid_ij(:,2));
ix_vec = gridix(grid_lin_idx); 

% Initialize output
ph_uw = nan(n_ps_grid, n_ifg, 'single');
valid_mask = (ix_vec > 0);

if any(valid_mask)
    % Extract valid components
    ph_smooth = uu.ph_uw(ix_vec(valid_mask), :); % Low-freq component
    ph_raw    = uw.ph_in(valid_mask, :);         % High-freq component
    
    % Combine: Unwrapped = Smooth_Grid + Angle(Residual)
    if isreal(uw.ph_in)
        ph_resid = exp(1i * (ph_raw - ph_smooth));
    else
        ph_resid = ph_raw .* exp(-1i * ph_smooth);
    end
    ph_uw(valid_mask, :) = ph_smooth + angle(ph_resid);
end

% -- Apply Predefined Phase Corrections (Cycle Skipping) --
if ~isempty(uw.ph_in_predef)
    predef_ix = ~isnan(uw.ph_in_predef);
    
    % Calculate global mean difference for alignment
    diff_val = ph_uw - uw.ph_in_predef;
    if exist('mean', 'builtin')
        meandiff = mean(diff_val, 'all', 'omitnan'); 
    else
        meandiff = nanmean(diff_val(:)); 
    end
    
    % Apply 2pi integer shift
    cycle_shift = 2 * pi * round(meandiff / (2*pi));
    uw.ph_in_predef = uw.ph_in_predef + cycle_shift;
    ph_uw(predef_ix) = uw.ph_in_predef(predef_ix);
end

msd = uu.msd;
fprintf('HPC Processing Complete. Total time: %.2f s\n', toc);

end