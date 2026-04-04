function [] = ps_unwrap()
%PS_UNWRAP Main driver for 3D Phase Unwrapping (HPC Block-IO Optimized).
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       2.0 (Block-IO Architecture)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Row-Block Streaming & IO Hand-off: Eradicated OOM crashes by streaming 
%      phase data in chunks to a temporary file ('ph_w_tmp.mat'), eliminating 
%      the need to pass massive matrices in memory.
%   2. Implicit Expansion: Replaced memory-heavy 'repmat' calls with modern 
%      implicit expansion, preventing massive temporary RAM spikes.
%   3. Structural Optimization: Refactored legacy code into clean blocks and 
%      merged 'sb_identify_good_pixels' to significantly reduce I/O overhead.
%   4. Dynamic Parameterization: Replaced hardcoded parameters with dynamic, 
%      platform-specific calculations for 'n_trial_wraps' (e.g., S1, ALOS, TSX).
%   5. Dependency Cleanup: Removed deprecated 'uw_nosnaphu' support, enforcing 
%      a robust Snaphu-based ('uw_3d') workflow for HPC consistency.
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
    small_baseline_flag = getparm('small_baseline_flag', 1);
    unwrap_patch_phase  = getparm('unwrap_patch_phase', 1);
    scla_deramp         = getparm('scla_deramp', 1);
    subtr_tropo         = getparm('subtr_tropo', 1);
    aps_name_param      = getparm('tropo_method', 1);
    drop_ifg_index      = getparm('drop_ifg_index', 1);
    
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

    ps = load(psname);
    unwrap_ifg_index = setdiff(1:ps.n_ifg, drop_ifg_index);

    % =====================================================================
    % 2. HPC Setup: Block Size & Temporary Storage
    % =====================================================================
    row_block_size = 1000000; % 1 Million points per block
    n_blocks = ceil(ps.n_ps / row_block_size);
    
    tmp_ph_file = 'ph_w_tmp.mat';
    m_tmp = matfile(tmp_ph_file, 'Writable', true);
    % Pre-allocate disk space quickly
    m_tmp.ph_w(ps.n_ps, ps.n_ifg) = complex(single(0), single(0));

    % Pre-allocate options.ph_uw_predef on disk if needed
    unwrap_hold_good = getparm('unwrap_hold_good_values', 1);
    if ~is_sbas || ~exist(phuwname, 'file')
        unwrap_hold_good = 'n'; 
    end
    if strcmpi(unwrap_hold_good, 'y')
        m_tmp.ph_uw_predef(ps.n_ps, ps.n_ifg) = single(NaN);
    end

    % =====================================================================
    % 3. Pre-load Global Variables (1D vectors or light matrices)
    % =====================================================================
    % BP
    if exist(['./', bpname, '.mat'], 'file')
        bp = load(bpname); bperp_mat = bp.bperp_mat;
    else
        bperp_vec = ps.bperp;
        if ~is_sbas, bperp_vec = bperp_vec([1:ps.master_ix-1, ps.master_ix+1:end]); end
        bperp_mat = repmat(bperp_vec', ps.n_ps, 1);
    end
    if ~is_sbas && size(bperp_mat, 2) ~= ps.n_ifg
        bperp_mat = [bperp_mat(:, 1:ps.master_ix-1), zeros(ps.n_ps, 1, 'single'), bperp_mat(:, ps.master_ix:end)];
    end

    % SCLA & APS (Flags and loading)
    scla_subtracted_sw = 0; ramp_subtracted_sw = 0;
    if exist([sclaname, '.mat'], 'file')
        scla = load(sclaname);
        if size(scla.K_ps_uw, 1) == ps.n_ps
            scla_subtracted_sw = 1;
            if strcmpi(scla_deramp, 'y') && isfield(scla, 'ph_ramp') && size(scla.ph_ramp, 1) == ps.n_ps
                ramp_subtracted_sw = 1;
            end
        else
            fprintf('   [Warning] SCLA point count mismatch. Skipping.\n');
            clear scla; scla_subtracted_sw = 0;
        end
    end

    aps_subtracted_sw = 0;
    if exist([apsname, '.mat'], 'file') && strcmpi(subtr_tropo, 'y')
        aps = load(apsname);
        [aps_corr, ~, ~] = ps_plot_tca(aps, aps_name_param); 
        aps_subtracted_sw = 1; clear aps;
    end

    % Phase Closure Pre-load
    if strcmpi(unwrap_hold_good, 'y')
        uw_prev_m = matfile(phuwname);
        loopname = ['phuw_loops', num2str(psver)];
        if ~exist([loopname,'.mat'], 'file'); sb_make_closure_loops; end
        load(loopname, 'intfg_loops'); 
        m_tmp.good_pixels(ps.n_ps, ps.n_ifg) = false;
    end

    % Data Source (Patch or RC)
    if strcmpi(unwrap_patch_phase, 'y')
        m_data = matfile(pmname); data_var = 'ph_patch';
    else
        m_data = matfile(rcname); data_var = 'ph_rc';
        if exist(['./', pmname, '.mat'], 'file')
            pm_topo = load(pmname, 'K_ps'); has_topo = isfield(pm_topo, 'K_ps');
        else
            has_topo = false;
        end
    end

    % =====================================================================
    % 4. BLOCK-IO Phase Preparation Loop (Subtracting Nuisances)
    % =====================================================================
    fprintf('   Pre-processing wrapped phase (Block-IO)...\n');
    for b = 1:n_blocks
        r_start = (b-1)*row_block_size + 1;
        r_end = min(b*row_block_size, ps.n_ps);
        rows = r_start:r_end;
        
        % A. Load Phase Chunk
        chunk = m_data.(data_var)(rows, :);
        
        % Align Master if needed
        if strcmpi(unwrap_patch_phase, 'y') && ~is_sbas
            chunk = [chunk(:, 1:ps.master_ix-1), ones(length(rows), 1, 'single'), chunk(:, ps.master_ix:end)];
        elseif ~strcmpi(unwrap_patch_phase, 'y') && has_topo
            % Add Topo back to RC
            % Implicit expansion used here: K_ps(rows) .* bperp_mat(rows,:)
            chunk = chunk .* exp(1i * (pm_topo.K_ps(rows) .* bperp_mat(rows, :)));
        end
        
        % B. Normalize Magnitude (Safe on chunk size)
        nz_mask = chunk ~= 0;
        chunk(nz_mask) = chunk(nz_mask) ./ abs(chunk(nz_mask));
        
        % C. Phase Closure Check (SBAS Hold Good)
        if strcmpi(unwrap_hold_good, 'y')
            uw_chunk = uw_prev_m.ph_uw(rows, :);
            closure_phasor = uw_chunk * intfg_loops'; 
            avg_resid = mean(closure_phasor, 1, 'omitnan'); 
            bias = 2 * pi * round(avg_resid / (2*pi));
            phase_resid = closure_phasor - bias;
            
            is_consistent = abs(phase_resid) <= 1;
            good_chunk = false(length(rows), ps.n_ifg);
            for m = 1:size(intfg_loops, 1)
                valid_idx = is_consistent(:, m);
                ifgs_in_loop = intfg_loops(m, :) ~= 0;
                good_chunk(valid_idx, ifgs_in_loop) = true;
            end
            
            predef_chunk = nan(length(rows), ps.n_ifg, 'single');
            predef_chunk(good_chunk) = uw_chunk(good_chunk);
            
            m_tmp.good_pixels(rows, :) = good_chunk;
        else
            predef_chunk = nan(length(rows), ps.n_ifg, 'single');
        end
        
        % D. Subtract SCLA and Ramps
        if scla_subtracted_sw
            % Implicit expansion: No repmat!
            chunk = chunk .* exp(-1i * (scla.K_ps_uw(rows) .* bperp_mat(rows, :)));
            if ~is_sbas
                chunk = chunk .* exp(-1i * scla.C_ps_uw(rows));
            end
            if ramp_subtracted_sw
                chunk = chunk .* exp(-1i * scla.ph_ramp(rows, :));
            end
            
            % Update Predefined Phase
            if strcmpi(unwrap_hold_good, 'y')
                predef_chunk = predef_chunk - (scla.K_ps_uw(rows) .* bperp_mat(rows, :));
                if ramp_subtracted_sw
                     predef_chunk = predef_chunk - scla.ph_ramp(rows, :);
                end
            end
        end
        
        % E. Subtract APS
        if aps_subtracted_sw
            chunk = chunk .* exp(-1i * aps_corr(rows, :));
            if strcmpi(unwrap_hold_good, 'y')
                predef_chunk = predef_chunk - aps_corr(rows, :);
            end
        end
        
        % F. Save Processed Chunk to Disk
        m_tmp.ph_w(rows, :) = chunk;
        if strcmpi(unwrap_hold_good, 'y')
            m_tmp.ph_uw_predef(rows, :) = predef_chunk;
        end
        
        if mod(b, 5) == 0 || b == n_blocks
            fprintf('      Progress: Block %d / %d pre-processed.\n', b, n_blocks);
        end
    end

    if strcmpi(unwrap_hold_good, 'y')
        good_pixels = m_tmp.good_pixels;
        save(goodname, 'good_pixels'); 
        clear good_pixels intfg_loops
    end

    % =====================================================================
    % 5. Unwrapping Options Configuration
    % =====================================================================
    options = struct('master_day', ps.master_day);
    options.time_win      = getparm('unwrap_time_win', 1);
    options.unwrap_method = getparm('unwrap_method', 1);
    options.grid_size     = getparm('unwrap_grid_size', 1);
    options.prefilt_win   = getparm('unwrap_gold_n_win', 1);
    options.goldfilt_flag = getparm('unwrap_prefilter_flag', 1);
    options.gold_alpha    = getparm('unwrap_gold_alpha', 1);
    options.la_flag       = getparm('unwrap_la_error_flag', 1);
    options.scf_flag      = getparm('unwrap_spatial_cost_func_flag', 1);

    % Calculate n_trial_wraps dynamically
    max_topo_err = getparm('max_topo_err', 1);
    lambda       = getparm('lambda', 1);
    platform_name = getparm('platform');
    if isempty(platform_name); platform_name = 'UNKNOWN'; end
    
    switch upper(platform_name)
        case {'SENTINEL-1', 'S1', 'S1A', 'S1B', 'SENTINEL'}
            default_rho = 880000; default_inc_rad = 37.5 * pi/180; la_offset = 4.5 * pi/180;
        case {'ALOS', 'ALOS2', 'ALOS-2'}
            default_rho = 760000; default_inc_rad = 34.3 * pi/180; la_offset = 4.5 * pi/180;
        otherwise
            default_rho = 630000; default_inc_rad = 35.0 * pi/180; la_offset = 5.0 * pi/180;
    end

    if isfield(ps, 'mean_incidence') && ~isempty(ps.mean_incidence)
        inc_mean = ps.mean_incidence;
    elseif exist(incname, 'file')
        inc = load(incname); inc_mean = mean(inc.inc(inc.inc ~= 0)); clear inc;
    else
        inc_mean = default_inc_rad;
    end

    rho = default_rho;
    if isfield(ps, 'mean_range') && ~isempty(ps.mean_range), rho = ps.mean_range; end

    max_K = max_topo_err / (lambda * rho * sin(inc_mean) / 4 / pi);
    options.n_trial_wraps = ((max(ps.bperp) - min(ps.bperp)) * max_K / (2 * pi));
    logit(sprintf('   Calculated n_trial_wraps = %f', options.n_trial_wraps));

    if is_sbas
        options.lowfilt_flag = 'n'; 
        ifgday_ix_use = ps.ifgday_ix;
        day_use = ps.day - ps.master_day;
    else
        ifgday_ix_use = [ones(ps.n_ifg, 1) * ps.master_ix, (1:ps.n_ifg)'];
        master_ix_in_day = sum(ps.master_day > ps.day) + 1;
        unwrap_ifg_index = setdiff(unwrap_ifg_index, master_ix_in_day); 
        day_use = ps.day - ps.master_day;
    end

    % Pass the Predefined array via matfile access directly inside uw_3d later
    % But we tell options we have it
    if strcmpi(unwrap_hold_good, 'y')
         options.ph_uw_predef_file = tmp_ph_file; % Pointer for uw_3d to load blocks
    end

    % =====================================================================
    % 6. EXECUTE UNWRAPPING (uw_3d)
    % =====================================================================
    % CRITICAL CHANGE: We pass the STRING 'ph_w_tmp.mat' instead of the matrix!
    [ph_uw_some, msd_some] = uw_3d(tmp_ph_file, ps.xy, day_use, ...
                                   ifgday_ix_use(unwrap_ifg_index, :), ...
                                   ps.bperp(unwrap_ifg_index), options, unwrap_ifg_index);

    % =====================================================================
    % 7. Post-processing: Add Back Nuisance Terms (Block-IO)
    % =====================================================================
    fprintf('   Post-processing and adding back nuisance terms (Block-IO)...\n');
    ph_uw = zeros(ps.n_ps, ps.n_ifg, 'single');
    
    if exist('msd_some', 'var')
        msd = zeros(ps.n_ifg, 1, 'single');
        msd(unwrap_ifg_index) = msd_some;
    end

    if strcmpi(unwrap_patch_phase, 'y')
        m_pm = matfile(pmname);
    end
    if exist([rcname, '.mat'], 'file')
        m_rc = matfile(rcname);
    end

    for b = 1:n_blocks
        r_start = (b-1)*row_block_size + 1;
        r_end = min(b*row_block_size, ps.n_ps);
        rows = r_start:r_end;
        
        uw_chunk = zeros(length(rows), ps.n_ifg, 'single');
        uw_chunk(:, unwrap_ifg_index) = ph_uw_some(rows, :);
        
        % Add back SCLA
        if scla_subtracted_sw
            uw_chunk = uw_chunk + (scla.K_ps_uw(rows) .* bperp_mat(rows, :)); 
            if ~is_sbas
                 uw_chunk = uw_chunk + scla.C_ps_uw(rows); 
            end
            if ramp_subtracted_sw
                uw_chunk = uw_chunk + scla.ph_ramp(rows, :); 
            end
        end
        
        % Add back APS
        if aps_subtracted_sw
            uw_chunk = uw_chunk + aps_corr(rows, :);
        end
        
        % Add back residual wrapped phase
        if strcmpi(unwrap_patch_phase, 'y')
            pm_patch_chunk = m_pm.ph_patch(rows, :);
            ph_w_patch = pm_patch_chunk ./ abs(pm_patch_chunk);
            if ~is_sbas
                ph_w_patch = [ph_w_patch(:, 1:ps.master_ix-1), zeros(length(rows), 1), ph_w_patch(:, ps.master_ix:end)];
            end
            
            rc_chunk = m_rc.ph_rc(rows, :);
            uw_chunk = uw_chunk + angle(rc_chunk .* conj(ph_w_patch));
        end
        
        uw_chunk(:, setdiff(1:ps.n_ifg, unwrap_ifg_index)) = 0;
        ph_uw(rows, :) = uw_chunk;
    end

    % =====================================================================
    % 8. Save Results & Cleanup
    % =====================================================================
    fprintf('   Saving results to %s...\n', phuwname);
    stamps_save(phuwname, ph_uw, msd);
    
    if exist(tmp_ph_file, 'file')
        delete(tmp_ph_file); % Clean up temporary file
    end
    
    logit(1); % Success

end