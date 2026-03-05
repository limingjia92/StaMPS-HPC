function [ph_all, ph_ramp] = ps_deramp(ps, ph_all)
%PS_DERAMP Removes orbital ramps from interferograms (HPC Optimized).
%   Estimates and subtracts orbital phase ramps and long-wavelength 
%   atmospheric signals using a robust, vectorized least-squares approach.
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
%   1. Vectorized Solver: Replaced iterative 'lscov' with vectorized QR decomposition 
%      (\), leveraging BLAS Level 3 for speedup and numerical stability.
%   2. Geometric Consistency: Implemented "Common Valid Pixel" strategy to ensure 
%      orbital planes are estimated from a consistent scatterer set across time.
%   3. Parameter Modernization: Replaced legacy 'load' with 'getparm' for 'deramp_degree'
%      and integrated spatial masking with coordinate rotation to exclude deformation.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, Septermber 2013
%   ======================================================================

    fprintf('PS_DERAMP: Estimating orbital ramps (HPC Optimized)...\n');
    
    % =====================================================================
    % 1. Parameter Retrieval
    % =====================================================================
    degree = getparm('deramp_degree');
    if isempty(degree), degree = 1; end
    
    % =====================================================================
    % 2. Spatial Masking & Exclusion Logic
    % =====================================================================
    mask_lon = getparm('deramp_mask_lon');
    mask_lat = getparm('deramp_mask_lat');
    
    % Initialize mask (false = keep, true = exclude)
    exclude_mask = false(ps.n_ps, 1); 
    
    if ~isempty(mask_lon) && ~isempty(mask_lat)
        % Ensure closed polygon
        if mask_lon(1) ~= mask_lon(end), mask_lon = [mask_lon mask_lon(1)]; end
        if mask_lat(1) ~= mask_lat(end), mask_lat = [mask_lat mask_lat(1)]; end
        
        % Convert mask to local coordinates (m)
        ll0 = ps.ll0;
        mask_xy_poly = llh2local([mask_lon' mask_lat']', ll0) * 1000;
        
        % Coordinate rotation to align with satellite track/heading
        heading = getparm('heading');
        if isempty(heading), heading = 0; end
        theta = (180 - heading) * pi / 180;
        rotm = [cos(theta), sin(theta); -sin(theta), cos(theta)];
        
        % Rotate mask polygon and PS points for check
        mask_xy_rot = (rotm * mask_xy_poly)';
        ps_xy_rot = (rotm * ps.xy(:, 2:3)')'; 
        
        % Identification of points inside the exclude polygon
        exclude_mask = inpolygon(ps_xy_rot(:,1), ps_xy_rot(:,2), mask_xy_rot(:,1), mask_xy_rot(:,2));
        
    end

    % =====================================================================
    % 3. Construct Design Matrix (A)
    % =====================================================================
    % Normalize coordinates to km to improve matrix conditioning
    x = ps.xy(:,2) / 1000; 
    y = ps.xy(:,3) / 1000; 
    ones_vec = ones(ps.n_ps, 1);
    

    % x axis means perpendicular to orbit direction, y axis means parallel to orbit direction.
    % so x axis is nearly east-west direction, y axis is nearly north-south direction 

    switch degree
        case 0,   A = [ones_vec];                                   % z = c
            fprintf('**** z = c\n');
        case 0.4, A = [x, ones_vec];                                % z = ax + c
            fprintf('**** z = ax+ c\n')
        case 0.6, A = [y, ones_vec];                                % z = by + c
            fprintf('**** z = by+ c\n')
        case 1,   A = [x, y, ones_vec];                             % z = ax + by + c
            fprintf('**** z = ax + by+ c\n')
        case 1.5, A = [x, y, x.*y, ones_vec];                       % z = ax + by + cxy + d
            fprintf('**** z = ax + by+ cxy + d\n')
        case 2,   A = [x.^2, y.^2, x.*y, ones_vec];                 % Simplified Quadratic
            fprintf('**** z = ax^2 + by^2 + cxy + d \n')
        case 2.5, A = [x.^2, y.^2, x.*y, x, y, ones_vec];           % Full Quadratic
            fprintf('**** z = ax^2 + by^2 + cxy + dx + ey +f \n')
        case 3,   A = [x.^3, y.^2.*x, y.^2.*x, x.^2, x.*y, y.^2, x.*y, ones_vec]; % Cubic (Simplified)
            fprintf('**** z = ax^3 + by^3 + cx^2y + dy^2x + ex^2 + fy^2 + gxy + h \n')
        otherwise, error('PS_DERAMP: Unsupported deramp degree.');
    end
    
    % =====================================================================
    % 4. Solver Strategy (Vectorized QR Decomposition)
    % =====================================================================
    
    % Identify "Common Valid Pixels" (Non-NaN in ALL ifgs AND Not Masked)
    nan_mask = isnan(ph_all);
    global_valid = ~any(nan_mask, 2) & ~exclude_mask; 
    
    n_valid = sum(global_valid);
    n_params = size(A, 2);
    n_ifg = size(ph_all, 2);

    % Threshold: Ensure enough points for a robust fit (e.g., > 2 * num_params)
    if n_valid > n_params * 2
        fprintf('   Solver: Using %d common valid points for vectorized estimation.\n', n_valid);
        
        % Global Solve: (N_valid x M_param) \ (N_valid x N_ifg)
        % Using MATLAB's mldivide (\) which employs QR decomposition for stability
        coeffs = A(global_valid, :) \ ph_all(global_valid, :); 
        
        % Reconstruct Ramp
        ph_ramp = A * coeffs;
        
        % Remove Ramp
        ph_all = ph_all - ph_ramp;
        
    else
        % Fallback Strategy: Loop through interferograms if NaN distribution varies significantly
        fprintf('   Solver: Common valid points insufficient. Falling back to per-interferogram loop.\n');
        
        ph_ramp = NaN(size(ph_all));
        
        for k = 1:n_ifg
            % Select valid points for current IFG
            valid_idx = ~nan_mask(:,k) & ~exclude_mask;
            
            if sum(valid_idx) > n_params + 5
                % Solve for current IFG
                c = A(valid_idx, :) \ ph_all(valid_idx, k);
                ramp_k = A * c;
                
                ph_ramp(:,k) = ramp_k;
                ph_all(:,k) = ph_all(:,k) - ramp_k;
            else
                 fprintf('   Warning: IFG %d skipped (insufficient valid points).\n', k);
            end
        end
    end
    
    fprintf('PS_DERAMP: Completed.\n');

end