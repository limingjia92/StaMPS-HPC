function [a, c1, c2, emp_vario] = ps_fit_vario(x, y, z, model, a0, c10, c20, decor_dist, Nlags, precalc_h_vec)
%PS_FIT_VARIO Estimate empirical variogram and fit theoretical model
%
%   Input:    - x               x-coordinates (or time array)
%             - y               y-coordinates (or zeros for time)
%             - z               z-values (phase residuals)
%             - model           Variogram model (2=exponential, 3=gaussian, 4=spherical)
%             - a0              Initial value range
%             - c10             Initial value sill
%             - c20             Initial value nugget
%             - decor_dist      Maximum distance to include
%             - Nlags           Number of lags for binning
%             - precalc_h_vec   (Optional) Pre-calculated distance vector
%
%   Output:   - a               Estimated range
%             - c1              Estimated sill
%             - c2              Estimated nugget
%             - emp_vario       Empirical variogram matrix
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
%   1. Vectorized Distances: Replaced double-for loops with native 'pdist' 
%      and added support for pre-calculated spatial distances to save memory.
%   2. Fast Binning: Replaced iterative 'find' calls with 'discretize' and 
%      'accumarray', reducing empirical variogram complexity to O(N).
%   3. Micro LM Solver: Completely bypassed 'lsqcurvefit' with a custom, 
%      inlined Levenberg-Marquardt solver using analytical Jacobians.
%   4. Robust Regularization: Added strict diagonal clamping and Tikhonov 
%      regularization to prevent singular matrix crashes on extreme data.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Freek van Leijen, November 2006
%   ======================================================================

% Ensure inputs are column vectors to maintain dimensional consistency for pdist
x = x(:);
y = y(:);
z = z(:);

c10 = double(c10);
c20 = double(c20);

% -------------------------------------------------------------------------
% 1. INITIALIZATION AND DISTANCE COMPUTATION
% -------------------------------------------------------------------------
% Utilize pre-calculated spatial distance vector if provided to completely
% eliminate redundant O(N^2) memory allocations across interferograms.
if nargin >= 10 && ~isempty(precalc_h_vec)
    h_vec = precalc_h_vec;
else
    coords = double([x, y]);
    h_vec = single(pdist(coords));
end

z_vec = single(pdist(double(z)));
gamma_vec = 0.5 * (z_vec.^2);

if ~exist('decor_dist', 'var') || isempty(decor_dist)
    max_dist = max(h_vec);
else
    max_dist = decor_dist;
end

if ~exist('Nlags', 'var') || isempty(Nlags)
    lags = 0:max_dist/50:max_dist;
    Nlags_orig = length(lags) - 1;
else
    lags = 0:max_dist/(Nlags-1):max_dist;
    Nlags_orig = Nlags - 1;
end

% -------------------------------------------------------------------------
% 2. EMPIRICAL VARIOGRAM ESTIMATION (FAST BINNING)
% -------------------------------------------------------------------------
% O(N) binning using discretize and accumarray, eliminating millions of 
% expensive find() loop iterations from the original implementation.
bins = discretize(h_vec, lags);
valid = ~isnan(bins);

if any(valid)
    bin_idx = bins(valid)';
    counts = accumarray(bin_idx, 1, [Nlags_orig, 1]);
    gamma_sums = accumarray(bin_idx, gamma_vec(valid)', [Nlags_orig, 1]);
    
    % Multiply counts by 2 to strictly match original i->j and j->i statistics
    counts = counts * 2; 
    gamma_mean = gamma_sums ./ max(counts/2, 1);
    
    lag_centers = (lags(1:end-1) + lags(2:end))' / 2;
    valid_bins = counts > 0;
    emp_vario = [counts(valid_bins), gamma_mean(valid_bins), lag_centers(valid_bins)];
else
    emp_vario = [];
end

% Truncate the empirical variogram to maintain consistency with the original heuristic
if size(emp_vario, 1) < round(Nlags_orig/1.75)
    last_lag = size(emp_vario, 1);
else
    last_lag = round(Nlags_orig/1.75);
end
emp_vario = emp_vario(1:last_lag, :);

if isempty(emp_vario)
    a = NaN; c1 = NaN; c2 = NaN;
    return;
end

if exist('decor_dist', 'var') && ~isempty(decor_dist)
    ind_max = find(emp_vario(:,3) < decor_dist, 1, 'last');
    if isempty(ind_max)
        ind_max = size(emp_vario, 1); 
    end
else
    ind_max = size(emp_vario, 1);
end

% -------------------------------------------------------------------------
% 3. THEORETICAL VARIOGRAM FITTING (CUSTOM LM SOLVER)
% -------------------------------------------------------------------------
H_obs = double(emp_vario(1:ind_max, 3));
y_obs = double(emp_vario(1:ind_max, 2));

% Initial guesses
a_est = double(a0);
c1_est = double(c10);
c2_est = double(c20);

max_iter = 50;
lambda = 0.01; % Levenberg-Marquardt damping factor

for iter = 1:max_iter
    % Step 3.1: Evaluate analytical model and Jacobian matrix
    if model == 2 % Exponential model
        y_est = c1_est * (1 - exp(-H_obs / a_est)) + c2_est;
        J_a = -c1_est * (H_obs / a_est^2) .* exp(-H_obs / a_est);
        J_c1 = 1 - exp(-H_obs / a_est);
        J_c2 = ones(size(H_obs));
        
    elseif model == 3 % Gaussian model
        y_est = c1_est * (1 - exp(-(H_obs / a_est).^2)) + c2_est;
        J_a = -2 * c1_est * (H_obs.^2 / a_est^3) .* exp(-(H_obs / a_est).^2);
        J_c1 = 1 - exp(-(H_obs / a_est).^2);
        J_c2 = ones(size(H_obs));
        
    elseif model == 4 % Spherical model
        h_a = min(H_obs / a_est, 1);
        y_est = c1_est * (1.5 * h_a - 0.5 * h_a.^3) + c2_est;
        J_a = zeros(size(H_obs));
        idx = H_obs < a_est;
        J_a(idx) = c1_est * (-1.5 * H_obs(idx) / a_est^2 + 1.5 * H_obs(idx).^3 / a_est^4);
        J_c1 = 1.5 * h_a - 0.5 * h_a.^3;
        J_c2 = ones(size(H_obs));
    else
        error('Unsupported variogram model in HPC version');
    end
    
    J = [J_a, J_c1, J_c2];
    residual = y_obs - y_est;
    
    % Step 3.2: Construct damped normal equations
    JTJ = J' * J;
    JTr = J' * residual;
    
    diag_JTJ = diag(JTJ);
    
    % Enforce strict regularization to avoid floating-point singularity
    diag_JTJ(diag_JTJ < 1e-10) = 1e-10; 
    
    % Levenberg-Marquardt damping with Tikhonov baseline
    H_lm = JTJ + lambda * diag(diag_JTJ) + 1e-8 * eye(3);
    
    % Suppress warnings for ill-conditioned matrices to prevent I/O lag in parfor
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:illConditionedMatrix');
    
    % Solve system using robust left division
    delta = H_lm \ JTr;
    
    warning('on', 'MATLAB:nearlySingularMatrix');
    warning('on', 'MATLAB:illConditionedMatrix');
    
    % Step 3.3: Update parameters with strict lower bounds
    a_new = max(a_est + delta(1), 1e-5);
    c1_new = max(c1_est + delta(2), 1e-5);
    c2_new = max(c2_est + delta(3), 1e-5);
    
    % Step 3.4: Evaluate step acceptance and adjust damping factor
    if model == 2
        y_new = c1_new * (1 - exp(-H_obs / a_new)) + c2_new;
    elseif model == 3
        y_new = c1_new * (1 - exp(-(H_obs / a_new).^2)) + c2_new;
    elseif model == 4
        h_a_new = min(H_obs / a_new, 1);
        y_new = c1_new * (1.5 * h_a_new - 0.5 * h_a_new.^3) + c2_new;
    end
    
    if sum((y_obs - y_new).^2) < sum(residual.^2)
        % Accept step, decrease damping (approach Gauss-Newton)
        a_est = a_new; c1_est = c1_new; c2_est = c2_new;
        lambda = max(lambda / 10, 1e-7);
        
        % Check for convergence
        if max(abs(delta ./ [a_est; c1_est; c2_est])) < 1e-4
            break;
        end
    else
        % Reject step, increase damping (approach Gradient Descent)
        lambda = min(lambda * 10, 1e7);
    end
end

% Finalize outputs
a = a_est;
c1 = c1_est;
c2 = c2_est;

end