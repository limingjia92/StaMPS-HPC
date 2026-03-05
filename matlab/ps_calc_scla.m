function [] = ps_calc_scla(use_small_baselines, coest_mean_vel)
%PS_CALC_SCLA Calculate SCLA and master atmosphere & orbit error (HPC Optimized)
%   Estimates spatially-correlated look angle error and master errors using 
%   vectorized GLS (L2) or IRLS (L1) approaches.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. L2 Vectorization & Precision: Replaced iterative 'lscov' with a vectorized GLS 
%      matrix operator, and fixed legacy bugs by strictly enforcing 'double' precision.
%   2. L1 Vectorized IRLS: Replaced the inefficient 'fminsearch' with fully vectorized 
%      IRLS, massively boosting execution speed and resolving local minima traps.
%   3. Architecture Refactoring: Restructured the core solver to strictly separate 
%      'scla_method' branches (L1/L2) and introduced block processing to prevent OOM.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, Nov 2006
%   ======================================================================

logit;
logit(sprintf('Estimating spatially-correlated look angle error (HPC Mode)...'), 2)

if nargin < 1, use_small_baselines = 0; end
if nargin < 2, coest_mean_vel = 0; end

% --- Parameter Loading ---
small_baseline_flag = getparm('small_baseline_flag', 1);
drop_ifg_index = getparm('drop_ifg_index', 1);
scla_method = getparm('scla_method', 1);
scla_deramp = getparm('scla_deramp', 1);
subtr_tropo = getparm('subtr_tropo', 1);
tropo_method = getparm('tropo_method', 1);

if isempty(scla_method), scla_method = 'L2'; end

if use_small_baselines ~= 0 && ~strcmpi(small_baseline_flag, 'y')
    error('   Use small baselines requested but there are none')
end

if use_small_baselines == 0
    scla_drop_index = getparm('scla_drop_index', 1);
else
    scla_drop_index = getparm('sb_scla_drop_index', 1);
    fprintf('   Using small baseline interferograms\n')
end

% --- File Paths ---
load psver
psname = ['./ps', num2str(psver)];
bpname = ['./bp', num2str(psver)];
meanvname = ['./mv', num2str(psver)];
ifgstdname = ['./ifgstd', num2str(psver)];
phuwsbresname = ['./phuw_sb_res', num2str(psver)];

if use_small_baselines == 0
    phuwname = ['./phuw', num2str(psver)];
    sclaname = ['./scla', num2str(psver)];
    apsname_old = ['./aps', num2str(psver)];
    apsname = ['./tca', num2str(psver)];
else
    phuwname = ['./phuw_sb', num2str(psver)];
    sclaname = ['./scla_sb', num2str(psver)];
    apsname_old = ['./aps_sb', num2str(psver)];
    apsname = ['./tca_sb', num2str(psver)];
end

if use_small_baselines == 0
    if exist([meanvname, '.mat'], 'file'), delete([meanvname, '.mat']); end
end

% --- Load Data ---
ps = load(psname);
if exist([bpname, '.mat'], 'file')
    bp = load(bpname);
else
    bperp = ps.bperp;
    if ~strcmpi(small_baseline_flag, 'y')
        bperp = bperp([1:ps.master_ix-1, ps.master_ix+1:end]);
    end
    bp.bperp_mat = repmat(bperp', ps.n_ps, 1);
end
uw = load(phuwname);

if strcmpi(small_baseline_flag, 'y') && use_small_baselines == 0
    unwrap_ifg_index = 1:ps.n_image;
else
    unwrap_ifg_index = setdiff(1:ps.n_ifg, drop_ifg_index);
end

% --- Preprocessing: Troposphere & Deramp ---
if strcmpi(subtr_tropo, 'y')
    if strcmpi(apsname, ['./tca', num2str(psver)])
        if strcmpi(getparm('small_baseline_flag'), 'y')
            sb_invert_tca(tropo_method);
        end
    end
    aps = load(apsname);
    [aps_corr, ~, tropo_method] = ps_plot_tca(aps, tropo_method);
    uw.ph_uw = uw.ph_uw - aps_corr;
end

ph_ramp = [];
if strcmpi(scla_deramp, 'y')
    fprintf('\n   deramping ifgs...\n')
    [~, ph_ramp] = ps_deramp(ps, uw.ph_uw);
    uw.ph_uw = uw.ph_uw - ph_ramp;
end

unwrap_ifg_index = setdiff(unwrap_ifg_index, scla_drop_index);

if exist([apsname_old, '.mat'], 'file')
    if strcmpi(subtr_tropo, 'y')
        fprintf('Warning: Removing atmosphere twice! Check subtr_tropo settings.\n')
    end
    aps = load(apsname_old);
    uw.ph_uw = uw.ph_uw - aps.ph_aps_slave;
end

% Reference Phase Correction
ref_ps = ps_setref(ps);
uw.ph_uw = uw.ph_uw - repmat(nanmean(uw.ph_uw(ref_ps, :), 1), ps.n_ps, 1);

% --- Prepare Bperp Matrix ---
if use_small_baselines == 0
    if strcmpi(small_baseline_flag, 'y')
        bperp_mat = zeros(ps.n_ps, ps.n_image, 'single');
        G = zeros(ps.n_ifg, ps.n_image);
        for i = 1:ps.n_ifg
            G(i, ps.ifgday_ix(i, 1)) = -1;
            G(i, ps.ifgday_ix(i, 2)) = 1;
        end
        if isfield(uw, 'unwrap_ifg_index_sm')
            unwrap_ifg_index = setdiff(uw.unwrap_ifg_index_sm, scla_drop_index);
        end
        unwrap_ifg_index = setdiff(unwrap_ifg_index, ps.master_ix);
        G = G(:, unwrap_ifg_index);
        bperp_mat(:, unwrap_ifg_index) = [G \ double(bp.bperp_mat')]';
    else
        bperp_mat = [bp.bperp_mat(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1, 'single'), bp.bperp_mat(:, ps.master_ix:end)];
    end
    day = diff(ps.day(unwrap_ifg_index));
    ph = double(diff(uw.ph_uw(:, unwrap_ifg_index), [], 2));
    bperp = diff(bperp_mat(:, unwrap_ifg_index), [], 2);
else
    bperp_mat = bp.bperp_mat;
    bperp = bperp_mat(:, unwrap_ifg_index);
    day = ps.ifgday(unwrap_ifg_index, 2) - ps.ifgday(unwrap_ifg_index, 1);
    ph = double(uw.ph_uw(:, unwrap_ifg_index));
end
clear bp

% --- Build Design Matrix (G) & Covariance ---
if coest_mean_vel == 0 || length(unwrap_ifg_index) < 4
    G = [ones(size(ph, 2), 1), double(mean(bperp)')];
else
    G = [ones(size(ph, 2), 1), double(mean(bperp)'), double(day)];
end

ifg_vcm = eye(ps.n_ifg);
if strcmpi(small_baseline_flag, 'y')
    if use_small_baselines == 0
        phuwres = load(phuwsbresname, 'sm_cov');
        if isfield(phuwres, 'sm_cov'), ifg_vcm = phuwres.sm_cov; end
    else
        phuwres = load(phuwsbresname, 'sb_cov');
        if isfield(phuwres, 'sb_cov'), ifg_vcm = phuwres.sb_cov; end
    end
else
    if exist([ifgstdname, '.mat'], 'file')
        ifgstd = load(ifgstdname);
        ifg_vcm = double(diag((ifgstd.ifg_std * pi / 180).^2));
    end
end

if use_small_baselines == 0
    ifg_vcm_use = eye(size(ph, 2));
else
    ifg_vcm_use = ifg_vcm(unwrap_ifg_index, unwrap_ifg_index);
end

while rcond(ifg_vcm_use) < 0.001
    ifg_vcm_use = ifg_vcm_use + eye(size(ifg_vcm_use)) * 0.01;
end

% =========================================================================
% CORE SOLVER: Estimate SCLA Parameters
% =========================================================================
block_size = 200000; % Block processing to prevent OOM
n_ps = ps.n_ps;
m_est = zeros(size(G, 2), n_ps);

fprintf('   Solving for SCLA using method: %s\n', scla_method);

switch upper(scla_method)
    case 'L2'
        % --- Method 1: L2 Norm (Vectorized GLS) ---
        W = inv(ifg_vcm_use);
        G_doub = double(G);
        H_L2 = (G_doub' * W * G_doub) \ (G_doub' * W);
        
        for b = 1:ceil(n_ps / block_size)
            idx = (b-1)*block_size + 1 : min(b*block_size, n_ps);
            d_chunk = ph(idx, :)'; 
            m_est(:, idx) = H_L2 * d_chunk;
        end
        
    case 'L1'
        % --- Method 2: L1 Norm (Vectorized IRLS with Damping) ---
        W = inv(ifg_vcm_use);
        G_doub = double(G);
        H_L2 = (G_doub' * W * G_doub) \ (G_doub' * W);
        
        % Column scaling for numerical stability
        scale_factors = max(abs(G_doub), [], 1);
        scale_factors(scale_factors == 0) = 1; 
        G_scaled = G_doub ./ scale_factors; 
        
        max_iter = 50;     
        epsilon = 1e-5;    
        tol = 1e-4;        
        num_params = size(G_scaled, 2); 
        
        for b = 1:ceil(n_ps / block_size)
            idx = (b-1)*block_size + 1 : min(b*block_size, n_ps);
            d_chunk = double(ph(idx, :)');
            
            % Strict L2 initialization
            m_init = H_L2 * d_chunk; 
            m_chunk_scaled = m_init .* scale_factors'; 
            
            % IRLS Loop
            for iter = 1:max_iter 
                m_old_scaled = m_chunk_scaled;
                resids = d_chunk - G_scaled * m_old_scaled;
                
                % Charbonnier smoothing
                weights = 1 ./ sqrt(resids.^2 + epsilon^2);
                m_new_scaled = zeros(size(m_chunk_scaled));
                
                % Fully unrolled Cramer's rule for extreme speed
                if num_params == 3
                    G1 = G_scaled(:,1); G2 = G_scaled(:,2); G3 = G_scaled(:,3);
                    w_G1 = weights .* G1; w_G2 = weights .* G2; w_G3 = weights .* G3;
                    A11 = sum(G1 .* w_G1, 1); A12 = sum(G1 .* w_G2, 1); A13 = sum(G1 .* w_G3, 1);
                    A22 = sum(G2 .* w_G2, 1); A23 = sum(G2 .* w_G3, 1); A33 = sum(G3 .* w_G3, 1);
                    b1 = sum(w_G1 .* d_chunk, 1); b2 = sum(w_G2 .* d_chunk, 1); b3 = sum(w_G3 .* d_chunk, 1);
                    
                    det_A = A11.*(A22.*A33 - A23.^2) - A12.*(A12.*A33 - A13.*A23) + A13.*(A12.*A23 - A13.*A22);
                    m_new_scaled(1,:) = ((A22.*A33 - A23.^2).*b1 - (A12.*A33 - A13.*A23).*b2 + (A12.*A23 - A13.*A22).*b3) ./ det_A;
                    m_new_scaled(2,:) = (-(A12.*A33 - A13.*A23).*b1 + (A11.*A33 - A13.^2).*b2 - (A11.*A23 - A12.*A13).*b3) ./ det_A;
                    m_new_scaled(3,:) = ((A12.*A23 - A13.*A22).*b1 - (A11.*A23 - A12.*A13).*b2 + (A11.*A22 - A12.^2).*b3) ./ det_A;
                elseif num_params == 2
                    G1 = G_scaled(:,1); G2 = G_scaled(:,2);
                    w_G1 = weights .* G1; w_G2 = weights .* G2;
                    A11 = sum(G1 .* w_G1, 1); A12 = sum(G1 .* w_G2, 1); A22 = sum(G2 .* w_G2, 1);
                    b1 = sum(w_G1 .* d_chunk, 1); b2 = sum(w_G2 .* d_chunk, 1);
                    det_A = A11.*A22 - A12.^2;
                    m_new_scaled(1,:) = (A22.*b1 - A12.*b2) ./ det_A;
                    m_new_scaled(2,:) = (-A12.*b1 + A11.*b2) ./ det_A;
                end
                
                % Damping update for stability
                m_chunk_scaled = 0.5 * m_old_scaled + 0.5 * m_new_scaled;
                
                if max(abs(m_chunk_scaled(:) - m_old_scaled(:))) < tol
                    break;
                end
            end
            
            m_est(:, idx) = m_chunk_scaled ./ scale_factors';
        end
        
    otherwise
        error('Unknown SCLA method. Use L2 or L1.');
end

K_ps_uw = m_est(2, :)';
if coest_mean_vel ~= 0 && size(m_est, 1) >= 3
    v_ps_uw = m_est(3, :)';
end

% --- Calculate Master Atmosphere & Orbit Error (C_ps_uw) ---
ph_scla = repmat(K_ps_uw, 1, size(bperp_mat, 2)) .* bperp_mat;

if use_small_baselines == 0
    unwrap_ifg_index = setdiff(unwrap_ifg_index, ps.master_ix);
    if coest_mean_vel == 0
        C_ps_uw = mean(uw.ph_uw(:, unwrap_ifg_index) - ph_scla(:, unwrap_ifg_index), 2);
    else
        G_master = [ones(length(unwrap_ifg_index), 1), ps.day(unwrap_ifg_index) - ps.day(ps.master_ix)];
        d_master = uw.ph_uw(:, unwrap_ifg_index) - ph_scla(:, unwrap_ifg_index);
        
        C_vcm = ifg_vcm(unwrap_ifg_index, unwrap_ifg_index);
        while rcond(C_vcm) < 0.001, C_vcm = C_vcm + eye(size(C_vcm))*0.01; end
        
        W_C = inv(C_vcm);
        G_m_doub = double(G_master);
        H_master = (G_m_doub' * W_C * G_m_doub) \ (G_m_doub' * W_C);
        
        m_master_est = zeros(2, n_ps);
        for b = 1:ceil(n_ps / block_size)
            idx = (b-1)*block_size + 1 : min(b*block_size, n_ps);
            d_chunk = double(d_master(idx, :)');
            m_master_est(:, idx) = H_master * d_chunk;
        end
        C_ps_uw = m_master_est(1, :)';
    end
else
    C_ps_uw = zeros(ps.n_ps, 1);
end

% --- Save Results ---
% HPC Opt: Removed legacy 'tmp_' backup logic to prevent disk bloat during iterative unwrapping
stamps_save(sclaname, ph_scla, K_ps_uw, C_ps_uw, ph_ramp, ifg_vcm)

logit(1);
end