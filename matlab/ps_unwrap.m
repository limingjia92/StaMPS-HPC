function [] = ps_unwrap()
%PS_UNWRAP Main driver for 3D Phase Unwrapping (HPC Optimized).
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
%   1. Structural Optimization: Refactored spaghetti code into clean, annotated 
%      blocks and merged 'sb_identify_good_pixels' to minimize I/O overhead.
%   2. Dynamic Parameterization: Replaced hardcoded values with platform-specific 
%      dynamic calculation for 'n_trial_wraps' (supporting S1, ALOS, TSX, etc.).
%   3. Dependency Cleanup: Removed deprecated 'uw_nosnaphu' support, enforcing 
%      'uw_3d' (Snaphu-based) workflow for consistency in HPC environments.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    logit;
    fprintf('Phase-unwrapping (HPC Version)...\n')

    % =====================================================================
    % 1. Configuration & File Setup
    % =====================================================================
    % Load parameters
    small_baseline_flag = getparm('small_baseline_flag', 1);
    unwrap_patch_phase  = getparm('unwrap_patch_phase', 1);
    scla_deramp         = getparm('scla_deramp', 1);
    subtr_tropo         = getparm('subtr_tropo', 1);
    aps_name_param      = getparm('tropo_method', 1);
    drop_ifg_index      = getparm('drop_ifg_index', 1);
    
    % Define file names based on processing version
    load psver
    psname  = ['ps', num2str(psver)];
    rcname  = ['rc', num2str(psver)];
    pmname  = ['pm', num2str(psver)];
    bpname  = ['bp', num2str(psver)];
    incname = ['inc', num2str(psver)];
    laname  = ['la', num2str(psver)];
    goodname= ['phuw_good', num2str(psver)];

    if strcmpi(small_baseline_flag, 'y')
        sclaname = ['scla_smooth_sb', num2str(psver)];
        apsname  = ['tca_sb', num2str(psver)];
        phuwname = ['phuw_sb', num2str(psver), '.mat'];
        is_sbas = true;
    else
        sclaname = ['scla_smooth', num2str(psver)];
        apsname  = ['tca', num2str(psver)];
        phuwname = ['phuw', num2str(psver), '.mat'];
        is_sbas = false;
    end

    % Load main PS structure
    ps = load(psname);
    unwrap_ifg_index = setdiff(1:ps.n_ifg, drop_ifg_index);

    % =====================================================================
    % 2. Baseline Matrix Preparation
    % =====================================================================
    if exist(['./', bpname, '.mat'], 'file')
        bp = load(bpname);
        bperp_mat = bp.bperp_mat;
    else
        % Construct on the fly if not saved
        bperp_vec = ps.bperp;
        if ~is_sbas
             % For PS, exclude master index from vector logic temporarily
             bperp_vec = bperp_vec([1:ps.master_ix-1, ps.master_ix+1:end]);
        end
        bperp_mat = repmat(bperp_vec', ps.n_ps, 1);
    end

    % Adjust Bperp matrix structure for PS (insert zero column for master)
    if ~is_sbas && size(bperp_mat, 2) ~= ps.n_ifg
        bperp_mat = [bperp_mat(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1, 'single'), bperp_mat(:, ps.master_ix:end)];
    end

    % =====================================================================
    % 3. Phase Data Loading (Patch vs RC)
    % =====================================================================
    if strcmpi(unwrap_patch_phase, 'y')
        fprintf('   Loading patch phase...\n');
        pm = load(pmname);
        ph_w = pm.ph_patch ./ abs(pm.ph_patch); % Normalize magnitude
        clear pm
        
        if ~is_sbas
            % Insert master column (zeros/ones) for PS
            ph_w = [ph_w(:, 1:ps.master_ix-1), ones(ps.n_ps, 1), ph_w(:, ps.master_ix:end)];
        end
    else
        fprintf('   Loading residual phase (RC)...\n');
        rc = load(rcname);
        ph_w = rc.ph_rc;
        clear rc;
        
        % Optionally remove topographic phase if PM exists
        if exist(['./', pmname, '.mat'], 'file')
            pm = load(pmname, 'K_ps');
            if isfield(pm, 'K_ps') && ~isempty(pm.K_ps)
                ph_w = ph_w .* exp(1i * (repmat(pm.K_ps, 1, ps.n_ifg) .* bperp_mat));
            end
        end
    end

    % Normalize phase (ensure unit magnitude)
    % Optimized: Matrix operation instead of double loop
    ph_w(ph_w ~= 0) = ph_w(ph_w ~= 0) ./ abs(ph_w(ph_w ~= 0));

    % =====================================================================
    % 4. Pre-processing: Subtract Nuisance Terms
    % =====================================================================
    options = struct('master_day', ps.master_day);
    scla_subtracted_sw = 0;
    ramp_subtracted_sw = 0;

    % 4.1 Recursive Unwrapping: Identify & Hold Good Pixels (SBAS Phase Closure)
    unwrap_hold_good = getparm('unwrap_hold_good_values', 1);
    
    % Only enable if SBAS mode and previous results exist
    if ~is_sbas || ~exist(phuwname, 'file')
        unwrap_hold_good = 'n'; 
    end
    
    if strcmpi(unwrap_hold_good, 'y')
        fprintf('   Identifying good pixels via Phase Closure check...\n');
        
        % 1. Load Data: Previous results and Loop definitions
        uw_prev = load(phuwname);
        loopname = ['phuw_loops', num2str(psver)];
        if ~exist([loopname,'.mat'], 'file'); sb_make_closure_loops; end
        load(loopname, 'intfg_loops'); 
        
        % 2. Calculate Closure Phase (Vectorized)
        % Project unwrapped phase onto closed loops via matrix multiplication
        closure_phasor = uw_prev.ph_uw * intfg_loops'; 
        
        % 3. Compute Residuals
        % Remove 2pi integer ambiguities to find the closure error
        avg_resid = mean(closure_phasor, 1, 'omitnan'); 
        bias = 2 * pi * round(avg_resid / (2*pi));
        phase_resid = closure_phasor - bias;
        
        % 4. Identify Consistent Loops (|Residual| <= 1 rad)
        is_consistent = abs(phase_resid) <= 1;
        
        % 5. Map Consistency to Interferograms
        % A pixel is 'good' for an IFG if valid in ANY loop containing that IFG
        good_pixels = false(size(uw_prev.ph_uw));
        for m = 1:size(intfg_loops, 1)
            valid_idx = is_consistent(:, m);
            ifgs_in_loop = intfg_loops(m, :) ~= 0;
            good_pixels(valid_idx, ifgs_in_loop) = true;
        end
        
        % 6. Set Predefined Phase Constraints
        options.ph_uw_predef = nan(size(ph_w), 'single');
        if ps.n_ps == size(good_pixels, 1) && ps.n_ps == size(uw_prev.ph_uw, 1)
            options.ph_uw_predef(good_pixels) = uw_prev.ph_uw(good_pixels);
            fprintf('   Phase closure check passed. Holding good pixels.\n');
        else
            fprintf('   [Warning] Mismatch in PS count. Skipping hold.\n');
        end
        
        save(goodname, 'good_pixels'); 
        clear uw_prev closure_phasor phase_resid good_pixels intfg_loops
    end

    % 4.2 Subtract SCLA (Spatially Correlated Look Angle Error) & Ramps
    if exist([sclaname, '.mat'], 'file')
        fprintf('   Subtracting SCLA and orbital ramps...\n');
        scla = load(sclaname);
        
        if size(scla.K_ps_uw, 1) == ps.n_ps
            scla_subtracted_sw = 1;
            
            % Subtract SCLA (K * Bperp)
            ph_w = ph_w .* exp(-1i * repmat(scla.K_ps_uw, 1, ps.n_ifg) .* bperp_mat);
            
            % For PS: Also subtract Master APS/AOE
            if ~is_sbas
                ph_w = ph_w .* repmat(exp(-1i * scla.C_ps_uw), 1, ps.n_ifg);
            end
            
            % Subtract Orbital Ramps
            if strcmpi(scla_deramp, 'y') && isfield(scla, 'ph_ramp') && size(scla.ph_ramp, 1) == ps.n_ps
                ramp_subtracted_sw = 1;
                ph_w = ph_w .* exp(-1i * scla.ph_ramp);
            end
            
            % Update predefined phase if holding values
            if strcmpi(unwrap_hold_good, 'y')
                options.ph_uw_predef = options.ph_uw_predef - repmat(scla.K_ps_uw, 1, ps.n_ifg) .* bperp_mat;
                if ramp_subtracted_sw
                     options.ph_uw_predef = options.ph_uw_predef - scla.ph_ramp;
                end
            end
        else
            fprintf('   [Warning] Wrong number of PS in SCLA file. Deleting and skipping.\n');
            delete([sclaname, '.mat']);
        end
        clear scla
    end
    
    % 4.3 Subtract Slave Atmosphere (APS)
    if exist([apsname, '.mat'], 'file') && strcmpi(subtr_tropo, 'y')
        fprintf('   Subtracting Slave APS...\n');
        aps = load(apsname);
        % Note: ps_plot_tca might be graphical, ensure it works in headless mode if needed
        [aps_corr, ~, ~] = ps_plot_tca(aps, aps_name_param); 
        
        ph_w = ph_w .* exp(-1i * aps_corr);
        
        if strcmpi(unwrap_hold_good, 'y')
            options.ph_uw_predef = options.ph_uw_predef - aps_corr;
        end
        clear aps
    end

    clear bp % Clean up large variables

    % =====================================================================
    % 5. Unwrapping Options Configuration
    % =====================================================================
    options.time_win      = getparm('unwrap_time_win', 1);
    options.unwrap_method = getparm('unwrap_method', 1);
    options.grid_size     = getparm('unwrap_grid_size', 1);
    options.prefilt_win   = getparm('unwrap_gold_n_win', 1);
    options.goldfilt_flag = getparm('unwrap_prefilter_flag', 1);
    options.gold_alpha    = getparm('unwrap_gold_alpha', 1);
    options.la_flag       = getparm('unwrap_la_error_flag', 1);
    options.scf_flag      = getparm('unwrap_spatial_cost_func_flag', 1);

    % --- Dynamic Parameter Calculation (Platform Specific) ---
    max_topo_err = getparm('max_topo_err', 1);
    lambda       = getparm('lambda', 1);
    platform_name = getparm('platform');
    
    if isempty(platform_name); platform_name = 'UNKNOWN'; end
    
    switch upper(platform_name)
        case {'SENTINEL-1', 'S1', 'S1A', 'S1B', 'SENTINEL'}
            default_rho = 880000; default_inc_rad = 37.5 * pi/180; la_offset = 4.5 * pi/180;
        case {'ALOS', 'ALOS2', 'ALOS-2'}
            default_rho = 760000; default_inc_rad = 34.3 * pi/180; la_offset = 4.5 * pi/180;
        case {'TERRASAR-X', 'TSX', 'TDX', 'PAZ'}
            default_rho = 630000; default_inc_rad = 35.0 * pi/180; la_offset = 5.0 * pi/180;
        case {'ENVISAT', 'ERS', 'ERS1', 'ERS2', 'ENV'}
            default_rho = 830000; default_inc_rad = 23.0 * pi/180; la_offset = 0.052;
        otherwise
            default_rho = 850000; default_inc_rad = 30.0 * pi/180; la_offset = 0.06;
    end

    % Determine Incidence Angle
    if isfield(ps, 'mean_incidence') && ~isempty(ps.mean_incidence)
        inc_mean = ps.mean_incidence;
    elseif exist(incname, 'file')
        inc = load(incname); inc_mean = mean(inc.inc(inc.inc ~= 0)); clear inc;
    elseif exist(laname, 'file')
        la = load(laname); inc_mean = mean(la(la ~= 0)) + la_offset; clear la;
    else
        inc_mean = default_inc_rad;
    end

    % Determine Range
    if isfield(ps, 'mean_range') && ~isempty(ps.mean_range)
        rho = ps.mean_range;
    else
        rho = default_rho;
    end

    % Calculate n_trial_wraps
    max_K = max_topo_err / (lambda * rho * sin(inc_mean) / 4 / pi);
    bperp_range = max(ps.bperp) - min(ps.bperp);
    options.n_trial_wraps = (bperp_range * max_K / (2 * pi));
    logit(sprintf('   Calculated n_trial_wraps = %f', options.n_trial_wraps));

    % Prepare vectors for uw_3d
    if is_sbas
        options.lowfilt_flag = 'n'; % SBAS usually disabled lowfilt
        ifgday_ix_use = ps.ifgday_ix;
        day_use = ps.day - ps.master_day;
    else
        % For PS, master is noise, remove it from list
        ifgday_ix_use = [ones(ps.n_ifg, 1) * ps.master_ix, (1:ps.n_ifg)'];
        master_ix_in_day = sum(ps.master_day > ps.day) + 1;
        unwrap_ifg_index = setdiff(unwrap_ifg_index, master_ix_in_day); 
        day_use = ps.day - ps.master_day;
    end

    if strcmpi(unwrap_hold_good, 'y')
        options.ph_uw_predef = options.ph_uw_predef(:, unwrap_ifg_index);
    end

    % =====================================================================
    % 6. EXECUTE UNWRAPPING (uw_3d)
    % =====================================================================
    % We force usage of uw_3d (HPC version uses Snaphu/3D logic)
    
    [ph_uw_some, msd_some] = uw_3d(ph_w(:, unwrap_ifg_index), ps.xy, day_use, ...
                                   ifgday_ix_use(unwrap_ifg_index, :), ...
                                   ps.bperp(unwrap_ifg_index), options);

    % Map results back to full size
    ph_uw = zeros(ps.n_ps, ps.n_ifg, 'single');
    msd = zeros(ps.n_ifg, 1, 'single');
    
    ph_uw(:, unwrap_ifg_index) = ph_uw_some;
    if exist('msd_some', 'var')
        msd(unwrap_ifg_index) = msd_some;
    end

    % =====================================================================
    % 7. Post-processing: Add Back Nuisance Terms
    % =====================================================================
    % Add back everything we subtracted in Step 4
    
    if scla_subtracted_sw
        fprintf('   Adding back SCLA and corrections...\n');
        scla = load(sclaname);
        
        ph_uw = ph_uw + (repmat(scla.K_ps_uw, 1, ps.n_ifg) .* bperp_mat); % Add SCLA
        
        if ~is_sbas
             ph_uw = ph_uw + repmat(scla.C_ps_uw, 1, ps.n_ifg); % Add Master APS
        end
        
        if ramp_subtracted_sw
            ph_uw = ph_uw + scla.ph_ramp; % Add Ramps
        end
        clear scla
    end

    if exist([apsname, '.mat'], 'file') && strcmpi(subtr_tropo, 'y')
        fprintf('   Adding back Slave APS...\n');
        aps = load(apsname);
        [aps_corr, ~, ~] = ps_plot_tca(aps, aps_name_param);
        ph_uw = ph_uw + aps_corr;
        clear aps
    end

    % If we unwrapped Patch Phase, add back Residual Phase
    if strcmpi(unwrap_patch_phase, 'y')
        pm = load(pmname);
        ph_w_patch = pm.ph_patch ./ abs(pm.ph_patch);
        clear pm
        
        if ~is_sbas
            ph_w_patch = [ph_w_patch(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1), ph_w_patch(:, ps.master_ix:end)];
        end
        
        rc = load(rcname);
        % This logic combines the unwrapped patch phase with the wrapped residual phase
        ph_uw = ph_uw + angle(rc.ph_rc .* conj(ph_w_patch));
    end

    % Ensure dropped interferograms are zero
    ph_uw(:, setdiff(1:ps.n_ifg, unwrap_ifg_index)) = 0;

    % =====================================================================
    % 8. Save Results
    % =====================================================================
    fprintf('   Saving results to %s...\n', phuwname);
    stamps_save(phuwname, ph_uw, msd);
    
    logit(1); % Success

end