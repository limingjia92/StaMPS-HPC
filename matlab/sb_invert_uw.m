function [] = sb_invert_uw()
%SB_INVERT_UW Invert unwrapped phase of short baseline ifgs 
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       2.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Two-Pass Block-IO Covariance: Eradicated >130GB OOM crashes by streaming 
%      massive RC and PM matrices through a 2-pass accumulator instead of full-RAM loading.
%   2. Streamed Matrix Inversion: Eliminated the 32GB full-matrix load of 'ph_uw_sb'. 
%      Phase data is now streamed chunk-by-chunk from disk directly into the solver.
%   3. Vectorized GLS Solver: Replaced the iterative 'lscov' pixel loop with a 
%      pre-calculated Generalized Least Squares (GLS) operator: 
%      H = (G'*inv(C)*G)^-1 * G'*inv(C).
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
    m_rc = matfile(rcname);
    m_pm = matfile(pmname);
    
    var_info = whos(m_pm);
    has_patch = any(strcmp({var_info.name}, 'ph_patch'));
    
    if has_patch
        fprintf('Calculating covariance matrix using Two-Pass Block-IO...\n');
        n_blocks = ceil(ps.n_ps / block_size);
        
        % Pass 1: Calculate Mean Noise per IFG
        sum_noise = zeros(1, ps.n_ifg, 'double');
        for b = 1:n_blocks
            idx_range = ((b - 1) * block_size + 1) : min(b * block_size, ps.n_ps);
            
            chunk_rc = m_rc.ph_rc(idx_range, :);
            chunk_pm = m_pm.ph_patch(idx_range, :);
            ph_noise = angle(chunk_rc .* conj(chunk_pm));
            
            sum_noise = sum_noise + double(sum(ph_noise, 1));
        end
        mean_noise = sum_noise / ps.n_ps;
        
        % Pass 2: Accumulate Covariance (X^T * X)
        cov_acc = zeros(ps.n_ifg, ps.n_ifg, 'double');
        for b = 1:n_blocks
            idx_range = ((b - 1) * block_size + 1) : min(b * block_size, ps.n_ps);
            
            chunk_rc = m_rc.ph_rc(idx_range, :);
            chunk_pm = m_pm.ph_patch(idx_range, :);
            ph_noise = double(angle(chunk_rc .* conj(chunk_pm)));
            
            % Mean Center the Chunk
            ph_noise = ph_noise - mean_noise;
            
            % Accumulate Cross-Products
            cov_acc = cov_acc + (ph_noise' * ph_noise);
        end
        
        sb_cov = cov_acc / (ps.n_ps - 1);
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

% =========================================================================
% CORE OPTIMIZATION: Generalized Least Squares (GLS) Operator
% =========================================================================
fprintf('Calculating inversion operator (GLS)...\n');

W = inv(C);
G2_doub = double(G2);
H_operator = (G2_doub' * W * G2_doub) \ (G2_doub' * W);

% --- Reference Phase Correction (Load only Ref Pixels) ---
m_phuwsb = matfile(phuwsbname);
ref_ps = ps_setref(ps);
ref_phase = nanmean(m_phuwsb.ph_uw(ref_ps, unwrap_ifg_index), 1);

% =========================================================================
% BLOCK PROCESSING (Streamed Matrix Inversion)
% =========================================================================
n_active_imgs = size(H_operator, 1);
ph_uw_inverted = zeros(ps.n_ps, n_active_imgs, 'single');

num_blocks = ceil(ps.n_ps / block_size);
fprintf('Starting streamed inversion on %d blocks...\n', num_blocks);

for b = 1:num_blocks
    start_idx = (b - 1) * block_size + 1;
    end_idx = min(b * block_size, ps.n_ps);
    idx_range = start_idx:end_idx;

    % 1. Stream Chunk Directly from Disk
    chunk_data = m_phuwsb.ph_uw(idx_range, unwrap_ifg_index);
    
    % 2. Apply Reference Phase Correction inline
    chunk_data = bsxfun(@minus, chunk_data, single(ref_phase));

    % 3. Apply Matrix Operator (High-performance BLAS level 3 operation)
    chunk_res = H_operator * double(chunk_data');

    % 4. Store Result
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