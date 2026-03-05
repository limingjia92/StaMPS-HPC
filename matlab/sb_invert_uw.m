function [] = sb_invert_uw()
%SB_INVERT_UW Invert unwrapped phase of short baseline ifgs 
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization Log:
%   1. Solver Replacement: Replaced iterative 'lscov' with pre-calculated
%      Generalized Least Squares (GLS) operator 'H'.
%      Math: m = (G'*inv(C)*G)^-1 * G'*inv(C) * d
%   2. Chunking: Implemented block processing (50k pixels/chunk) to minimize
%      RAM footprint during transpose operations 
%
% ======================================================================
%   ORIGINAL HEADER (StaMPS)
% ======================================================================
%   Original Author: Andy Hooper, September 2006
% ======================================================================


% --- Configuration ---
block_size = 200000; % Chunk size to optimize RAM/Cache usage

% --- Load Paths and Metadata ---
load psver
psname = ['./ps', num2str(psver)];
rcname = ['./rc', num2str(psver)];
pmname = ['./pm', num2str(psver)];
phuwsbname = ['./phuw_sb', num2str(psver)];
phuwsbresname = ['./phuw_sb_res', num2str(psver)];
phuwname = ['./phuw', num2str(psver)];

ps = load(psname);

% --- Covariance Matrix Setup ---
drop_ifg_index = getparm('drop_ifg_index');
unwrap_ifg_index = setdiff(1:ps.n_ifg, drop_ifg_index);

if exist([pmname, '.mat'], 'file')
    rc = load(rcname);
    pm = load(pmname);
    if isfield(pm, 'ph_patch') && ~isempty(pm.ph_patch)
        ph_noise = angle(rc.ph_rc .* conj(pm.ph_patch));
        sb_cov = double(cov(ph_noise));
        clear ph_noise pm rc
    else
        sb_cov = eye(ps.n_ifg);
    end
else
    sb_cov = eye(ps.n_ifg);
end

C = sb_cov(unwrap_ifg_index, unwrap_ifg_index);

% Regularization to ensure positive definiteness
while rcond(C) < 0.001
    C = C + eye(size(C, 1)) * 0.01;
end

% --- Design Matrix (G) Construction ---
G = zeros(ps.n_ifg, ps.n_image);
for i = 1:ps.n_ifg
    G(i, ps.ifgday_ix(i, 1)) = -1;
    G(i, ps.ifgday_ix(i, 2)) = 1;
end
G(:, ps.master_ix) = 0; % Remove master reference

G2 = G(unwrap_ifg_index, :);
nzc_ix = sum(abs(G2)) ~= 0;
G2 = G2(:, nzc_ix);

if rank(G2) < size(G2, 2)
    stamps_save(phuwsbresname, sb_cov);
    error('Error: Isolated subsets detected (network disconnected).');
end

% --- Load Phase Data ---
phuwsb = load(phuwsbname);
ph_uw_sb = phuwsb.ph_uw(:, unwrap_ifg_index);
clear phuwsb

% --- Reference Phase Correction ---
ref_ps = ps_setref(ps);
ref_phase = nanmean(ph_uw_sb(ref_ps, :), 1);
ph_uw_sb = bsxfun(@minus, ph_uw_sb, ref_phase);

% =========================================================================
% CORE OPTIMIZATION: Generalized Least Squares (GLS) Operator
% =========================================================================
fprintf('Calculating inversion operator (GLS)...\n');

% Original Logic (Legacy):
%   ph_uw = lscov(G2, ph_uw_sb', C);
%   This iteratively solved d = G*m for every single pixel, which is slow 
%   and memory-intensive for millions of points.

% HPC Optimization (Matrix Operator):
%   We solve the linear system d = G*m using the analytical GLS solution:
%       m_hat = (G' * C^-1 * G)^-1 * G' * C^-1 * d
%
%   Since G (Design Matrix) and C (Covariance) are constant for all pixels,
%   we pre-calculate the projection matrix 'H_operator':
%       H_operator = (G' * C^-1 * G)^-1 * G' * C^-1
%
%   The inversion then simplifies to a single matrix multiplication:
%       m_hat = H_operator * d
%   This leverages BLAS/LAPACK optimizations for massive speedup.

W = inv(C);
G2_doub = double(G2);
H_operator = (G2_doub' * W * G2_doub) \ (G2_doub' * W);

% =========================================================================
% BLOCK PROCESSING
% =========================================================================
% Process data in chunks to minimize RAM usage during transpose operations.
% Implicit multithreading in matrix multiplication replaces explicit parfor.

[n_ps, ~] = size(ph_uw_sb);
n_active_imgs = size(H_operator, 1);
ph_uw_inverted = zeros(n_ps, n_active_imgs, 'single');

num_blocks = ceil(n_ps / block_size);
fprintf('Starting inversion on %d blocks...\n', num_blocks);

for b = 1:num_blocks
    start_idx = (b - 1) * block_size + 1;
    end_idx = min(b * block_size, n_ps);
    idx_range = start_idx:end_idx;

    % 1. Extract Chunk (Transpose only the small block to double)
    chunk_data = double(ph_uw_sb(idx_range, :)');

    % 2. Apply Matrix Operator (High-performance BLAS level 3 operation)
    chunk_res = H_operator * chunk_data;

    % 3. Store Result
    ph_uw_inverted(idx_range, :) = single(chunk_res');
end

% --- Reconstruct & Save ---
ph_uw = zeros(ps.n_ps, ps.n_image, 'single');
ph_uw(:, nzc_ix) = ph_uw_inverted;
clear ph_uw_inverted

ph_res = single(G * ph_uw')';

sm_cov = zeros(ps.n_image);
sm_cov(nzc_ix, nzc_ix) = inv(G2_doub' * W * G2_doub);

unwrap_ifg_index_sm = 1:ps.n_image;
nzc_ix(ps.master_ix) = 1;
unwrap_ifg_index_sm = unwrap_ifg_index_sm(nzc_ix);

stamps_save(phuwname, ph_uw, unwrap_ifg_index_sm);
stamps_save(phuwsbresname, ph_res, sb_cov, sm_cov);

fprintf('SBAS Inversion Complete.\n');

logit(1);
end