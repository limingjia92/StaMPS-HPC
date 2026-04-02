function [] = ps_calc_ifg_std
% PS_CALC_IFG_STD (HPC Row-Block Optimized Version)
%   Calculate standard deviation (noise level) for each interferogram.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       2.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   OPTIMIZATION HIGHLIGHTS:
%   1. Mathematical Simplification: Replaced computationally expensive 
%      Complex Exponential operations (`exp(-j*...)`) with direct Real 
%      Number arithmetic (`phase - correction`).
%   2. Dimensional Inversion: Defeats HDF5 decompression bottlenecks by reading
%      data row-by-row (1 million points at a time across ALL interferograms) 
%      instead of column-by-column. Ensures single-pass sequential disk read.
%   3. Minimal RAM: Maximum memory strictly capped below ~10GB.
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

fprintf('\nEstimating noise standard deviation (degrees) using Row-Block IO...\n')

small_baseline_flag = getparm('small_baseline_flag');

load psver
psname = ['ps', num2str(psver)];
phname = ['ph', num2str(psver)];
pmname = ['pm', num2str(psver)];
bpname = ['bp', num2str(psver)];
ifgstdname = ['ifgstd', num2str(psver)];

if strcmpi(small_baseline_flag, 'y')
    ps = load(psname, 'xy', 'master_day', 'day', 'n_ifg', 'ifgday');
else
    ps = load(psname, 'xy', 'master_day', 'day', 'n_ifg');
end
n_ps = length(ps.xy);
master_ix = sum(ps.master_day > ps.day) + 1;
n_ifg = ps.n_ifg;

% -------------------------------------------------------------------------
% 1. Setup Memory-Mapped matfile Access
% -------------------------------------------------------------------------
if exist([phname, '.mat'], 'file')
    m_ph = matfile(phname);
else
    m_ph = matfile(psname); % Fallback
end
m_bp = matfile(bpname);
m_pm = matfile(pmname);

% Load 1D Vectors globally
K_ps = m_pm.K_ps;
if ~strcmpi(small_baseline_flag, 'y')
    C_ps = m_pm.C_ps;
end

% Check BP matrix shape
sz_bp = size(m_bp, 'bperp_mat');
is_bp_vector = (sz_bp(1) == 1);
if is_bp_vector
    bp_vec = m_bp.bperp_mat;
end

% Pre-allocate the Sum of Squares accumulator
sum_sq = zeros(1, n_ifg, 'single');

% --- Row-Block Setup ---
row_block_size = 1000000; % Process 1 million points per block
n_blocks = ceil(n_ps / row_block_size);

fprintf('   Processing %d points in %d Blocks (Block Size: %d rows)...\n', n_ps, n_blocks, row_block_size);

% -------------------------------------------------------------------------
% 2. Row-Block Processing Loop
% -------------------------------------------------------------------------
for b = 1:n_blocks
    row_start = (b-1)*row_block_size + 1;
    row_end = min(b*row_block_size, n_ps);
    rows = row_start:row_end;
    n_rows_current = length(rows);
    
    % 1) Fast Sequential Load: Read a thick horizontal slice of the matrix
    ph_angle_chunk = angle(m_ph.ph(rows, :));
    
    % 2) Extract PM & BP data and compute Phase Diff
    if strcmpi(small_baseline_flag, 'y')
        pm_angle_chunk = angle(m_pm.ph_patch(rows, :));
        
        if is_bp_vector
            bp_chunk = bp_vec; % (1 x n_ifg)
        else
            bp_chunk = m_bp.bperp_mat(rows, :);
        end
        
        topo_phase = K_ps(rows) .* bp_chunk; 
        ph_diff = wrap_phase(ph_angle_chunk - pm_angle_chunk - topo_phase);
        
    else
        % Non-SBAS logic: Align matrices by inserting master_ix dynamically
        pm_raw = m_pm.ph_patch(rows, :);
        pm_angle_chunk = [angle(pm_raw(:, 1:master_ix-1)), ...
                          zeros(n_rows_current, 1, 'single'), ...
                          angle(pm_raw(:, master_ix:end))];
        clear pm_raw;
        
        if is_bp_vector
            bp_chunk_full = [repmat(bp_vec(1, 1:master_ix-1), n_rows_current, 1), ...
                             zeros(n_rows_current, 1, 'single'), ...
                             repmat(bp_vec(1, master_ix:end), n_rows_current, 1)];
        else
            bp_raw = m_bp.bperp_mat(rows, :);
            bp_chunk_full = [bp_raw(:, 1:master_ix-1), ...
                             zeros(n_rows_current, 1, 'single'), ...
                             bp_raw(:, master_ix:end)];
            clear bp_raw;
        end
        
        term_correction = (K_ps(rows) .* bp_chunk_full) + C_ps(rows); 
        ph_diff = wrap_phase(ph_angle_chunk - pm_angle_chunk - term_correction);
    end
    
    % 3) Accumulate Sum of Squares for this row block
    sum_sq = sum_sq + sum(ph_diff.^2, 1);
    
    fprintf('      Progress: Block %d / %d processed.\n', b, n_blocks);
end

% Calculate final standard deviation
ifg_std = sqrt(sum_sq / n_ps) * 180 / pi;

ifg_std = ifg_std';

% -------------------------------------------------------------------------
% 3. Display Results and Save
% -------------------------------------------------------------------------
fprintf('\n');
if strcmpi(small_baseline_flag, 'y')
    for i = 1:n_ifg
        fprintf('%3d %s -- %s %3.2f\n', i, datestr(ps.ifgday(i,1)), datestr(ps.ifgday(i,2)), ifg_std(i));
    end
else
    for i = 1:n_ifg
        fprintf('%3d %s %3.2f\n', i, datestr(ps.day(i)), ifg_std(i));
    end
end
fprintf('\n')

stamps_save(ifgstdname, ifg_std); 

end

% --- Helper Function ---
function ph_wrapped = wrap_phase(ph)
    % Wrap phase to (-pi, pi]
    ph_wrapped = mod(ph + pi, 2*pi) - pi;
end