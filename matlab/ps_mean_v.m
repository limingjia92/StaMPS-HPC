function mean_v_std = ps_mean_v(ph_uw_mm, G, ifg_cov, n_boot)
%PS_MEAN_V Calculate standard deviations of mean velocities (HPC Optimized).
%   Takes fully assembled phase and design matrices to perform ultra-fast, 
%   Cholesky-prewhitened bootstrapping without any disk I/O.
%
%   INPUTS:
%   ph_uw_mm : [N_ifg x N_ps] Unwrapped phase converted to mm (Ready for LSQ)
%   G        : [N_ifg x 2] Design matrix [ones, delta_time_years]
%   ifg_cov  : [N_ifg x N_ifg] Covariance matrix (sm_cov or sb_cov)
%   n_boot   : Integer, number of bootstrap iterations (e.g., 1500)
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization & Algorithm Enhancements:
%   1. Architectural Decoupling: Bypassed redundant file I/O by directly 
%      accepting pre-assembled phase matrices from the dynamic assembler.
%   2. Mathematical Rigor: Applied Cholesky pre-whitening to enforce i.i.d. 
%      samples, fixing the severe statistical artifact in original bootstrapping.
%   3. Computational Speedup: Fed pre-whitened matrices into dual-parameter 
%      'lscov' (LAPACK QR) to eliminate massive memory reallocation overheads.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, Jan 2008
%   ======================================================================

    fprintf('Calculating standard deviation of mean velocity (HPC Optimized)...\n');
    
    n_ps = size(ph_uw_mm, 2);
    N_ifg = size(ph_uw_mm, 1);
    mean_v_std = NaN(n_ps, 1, 'single');
    
    % Apply Cholesky pre-whitening to convert Weighted Least Squares (WLS) 
    % to Ordinary Least Squares (OLS) for valid and fast bootstrapping.
    try
        L = chol(inv(ifg_cov)); 
        ph_white = L * ph_uw_mm;
        G_white = L * G;
    catch
        warning('Covariance matrix is singular, falling back to OLS without pre-whitening');
        ph_white = ph_uw_mm;
        G_white = G;
    end
    
    % Initialize random number generator
    rand('twister', sum(100*clock)); 
    
    % Process in chunks to prevent memory overflow
    chunk_size = 10000; 
    
    print_step = max(chunk_size, ceil((n_ps * 0.05) / chunk_size) * chunk_size);
    
    i = 1;
    while i <= n_ps
        i_end = min(i + chunk_size - 1, n_ps);
        ph_bit = ph_white(:, i:i_end);
        
        % Identify valid PS points (without NaN values)
        ix_non_nan = sum(isnan(ph_bit), 1) == 0;
        
        if sum(ix_non_nan) > 0
            if sum(ix_non_nan) == length(ix_non_nan)
                % Process chunk where all points are valid
                [mean_v_dist, ~] = bootstrp(n_boot, @(x) single(lscov(G_white(x,:), ph_bit(x,:))), (1:N_ifg)');                
                temp = std(mean_v_dist(:, 2:2:end))';
            else
                % Process chunk containing NaN values
                [mean_v_dist_temp, ~] = bootstrp(n_boot, @(x) single(lscov(G_white(x,:), ph_bit(x, ix_non_nan))), (1:N_ifg)');                
                temp = NaN(size(ph_bit, 2), 1); 
                temp(ix_non_nan) = std(mean_v_dist_temp(:, 2:2:end))';
            end
            mean_v_std(i:i_end) = temp;       
        end
        
        % 动态进度输出
        if mod(i_end, print_step) == 0 || i_end == n_ps
            fprintf('  -> %d / %d PS processed\n', i_end, n_ps);
        end
        
        i = i + chunk_size;
    end
end