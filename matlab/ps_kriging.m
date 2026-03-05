function [sol, var_sol, min_dist, mean_dist, nof_obs] = ps_kriging(x, x0, Nmax, kriging_method, cv_model)
%PS_KRIGING Perform spatial or temporal Kriging interpolation
%
%   Input:    - x               Data points [cn1 cn2 z] (subsetted neighbors)
%             - x0              Target point [cn1 cn2] (1x2 vector)
%             - Nmax            Maximum number of nearest neighbors to use
%             - kriging_method  1 = Simple Kriging, 2 = Ordinary Kriging
%             - cv_model        Covariance model matrix
%
%   Output:   - sol             Interpolated solution
%             - var_sol         Kriging variance
%             - min_dist        Minimum distance to a neighbor
%             - mean_dist       Mean distance to neighbors
%             - nof_obs         Number of observations within range
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
%   1. Data Truncation: Enforced Nmax truncation with rapid squared-distance 
%      sorting to prevent memory spikes during temporal Kriging.
%   2. Implicit Expansion: Eliminated memory-heavy 'repmat' calls using 
%      vectorized distance matrix calculations for superior throughput.
%   3. Inlined Covariance: Hardcoded variogram models directly into the 
%      loop, stripping millions of costly sub-function evaluation overheads.
%   4. Robust Solver: Utilized robust left-division (\) for symmetric 
%      systems and suppressed ill-conditioned warnings to prevent I/O lag.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Freek van Leijen, November 2006
%   ======================================================================

% -------------------------------------------------------------------------
% 1. DATA PREPARATION AND NEIGHBORHOOD TRUNCATION
% -------------------------------------------------------------------------
% Ensure the number of provided observations does not exceed Nmax.
% This acts as a strict memory safeguard, primarily for temporal Kriging 
% where points within a radius might be passed without prior truncation.
Nx = size(x, 1);
Nmax = min(Nmax, Nx);

if Nx > Nmax
    dx_sq = (x(:,1) - x0(1)).^2 + (x(:,2) - x0(2)).^2;
    [~, sort_idx] = sort(dx_sq);
    x = x(sort_idx(1:Nmax), :);
    Nx = Nmax;
end

% Ordinary Kriging requires at least 3 valid points to formulate constraints
if Nx < 3
    kriging_method = 2; 
end

% -------------------------------------------------------------------------
% 2. DISTANCE MATRIX COMPUTATION
% -------------------------------------------------------------------------
% Utilize implicit expansion to compute pairwise distances efficiently
% without the severe memory and time overhead associated with repmat.
X_coords = x(:, 1:2);

% Calculate distances from the target point to all neighboring data points
h0 = sqrt((X_coords(:,1) - x0(1)).^2 + (X_coords(:,2) - x0(2)).^2);

min_dist = min(h0);
mean_dist = mean(h0);
nof_obs = sum(h0 < cv_model(2,2)); 

% Calculate the pairwise distance matrix between all neighboring points
h = sqrt((X_coords(:,1) - X_coords(:,1)').^2 + (X_coords(:,2) - X_coords(:,2)').^2);

% -------------------------------------------------------------------------
% 3. COVARIANCE MATRIX ASSEMBLY
% -------------------------------------------------------------------------
% Inline the variogram models to evaluate covariance directly. This bypasses
% the overhead of calling external functions millions of times in parfor.
k = zeros(Nx, Nx, 'double');
k0 = zeros(Nx, 1, 'double');
Nmodel = size(cv_model, 1);

for w = 1:Nmodel
    model = cv_model(w, 1);
    a = cv_model(w, 2);
    c = cv_model(w, 3);
    
    if model == 1 % Nugget
        k = k + eye(Nx) * c;
    elseif model == 2 % Exponential
        cv = c * exp(-h / a);
        k = k + cv;
        k0 = k0 + c * exp(-h0 / a);
    elseif model == 3 % Gaussian
        cv = c * exp(-(h / a).^2);
        k = k + cv;
        k0 = k0 + c * exp(-(h0 / a).^2);
    elseif model == 4 % Spherical
        h_a = min(h / a, 1);
        cv = c * (1 - 1.5 * h_a + 0.5 * h_a.^3);
        k = k + cv;
        
        h0_a = min(h0 / a, 1);
        k0 = k0 + c * (1 - 1.5 * h0_a + 0.5 * h0_a.^3);
    end
end

% -------------------------------------------------------------------------
% 4. KRIGING SYSTEM SOLUTION
% -------------------------------------------------------------------------
% Mean center the observations to stabilize the numerical solution
avg = mean(x(:,3));
z = x(:,3) - avg;

% Temporarily disable nearly singular matrix warnings. For massive parfor 
% execution, console warning spam causes severe communication bottlenecks.
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

if kriging_method == 1 % Simple Kriging
    % Solve the symmetric positive definite system
    l = k \ k0;
    
    sol = l' * z + avg;
    var_sol = k(1,1) - l' * k0;
    
elseif kriging_method == 2 % Ordinary Kriging
    % Append the non-bias constraint row and column
    k_sys = [k, ones(Nx, 1); ones(1, Nx), 0];
    k0_sys = [k0; 1];
    
    % Solve the augmented system using robust left division (\)
    l = k_sys \ k0_sys;
    
    % Extract weights and compute the final interpolated solution
    weights = l(1:Nx);
    sol = weights' * z + avg;
    var_sol = k(1,1) - l' * k0_sys;
    
else
    % Fallback for unimplemented models
    sol = NaN; 
    var_sol = NaN;
    warning('HPC Kriging currently only supports Simple (1) or Ordinary (2) Kriging.');
end

% Restore warning states to maintain strict code execution standards
warning('on', 'MATLAB:nearlySingularMatrix');
warning('on', 'MATLAB:illConditionedMatrix');

if isnan(sol)
    var_sol = NaN;
    min_dist = NaN;
    mean_dist = NaN;
end

end