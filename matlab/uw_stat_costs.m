function [] = uw_stat_costs(options)
%UW_STAT_COSTS Generate statistical costs for MAP phase unwrapping (HPC Optimized).
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   =====================================================================
%   Author:        Mingjia Li
%   Date:          August 17, 2026
%   Version:       1.1
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Parallel Execution: Uses parfor to process interferograms concurrently
%      (one IFG per MATLAB worker).
%   2. SNAPHU Tiling: Uses a 3x3 tiled SNAPHU configuration with NPROC=1.
%      IFG-level parallelism is handled by MATLAB parfor, avoiding nested
%      SNAPHU process oversubscription.
%   3. SNAPHU Re-optimization: Retains the '-S' option. '-S' does NOT enable
%      tile mode; it performs a full single-tile re-optimization after the
%      multi-tile initialization/assembly. Tile mode itself is controlled by
%      NTILEROW/NTILECOL. A future version may expose '-S' through options.
%   4. Logic Optimization: Vectorized cost calculation while preserving the
%      StaMPS statistical-cost formulation.
%
%   Correctness & Robustness Fixes (v1.1):
%   5. SNAPHU Detection: Uses 'command -v snaphu' to test availability and
%      parses 'snaphu --info' only for version text, avoiding reliance on the
%      non-standard help/info exit status used by SNAPHU.
%   6. Version Policy: '-S' requires SNAPHU >= 2.0.0. Versions 2.0.0-2.0.5
%      are allowed with a strong production warning; v2.0.6 is allowed but
%      v2.0.7+ is recommended for the current tiled production workflow.
%   7. StaMPS Logic: Restores loading of predef_ix from uw_space_time.mat so
%      predefined low-cost edges are honored when present.
%   8. I/O Safety: Uses native byte order, validates file creation and SNAPHU
%      output dimensions, preserves failed job directories, and prevents a
%      partially failed unwrap from being silently saved as a valid product.
%   9. Pool Lifecycle: Only closes a parallel pool created by this function.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, May 2007
%   ======================================================================

if nargin < 1 || isempty(options)
    options = struct();
elseif ~isstruct(options)
    error('StaMPS_HPC:InvalidOptions', 'options must be a struct.');
end

if isfield(options, 'unwrap_method')
    unwrap_method = options.unwrap_method;
else
    unwrap_method = '3D';
end

if isfield(options, 'variance')
    variance = options.variance;
else
    variance = [];
end

% =========================================================================
% SNAPHU DETECTION AND VERSION POLICY
% =========================================================================
% Do not use the exit status from 'snaphu -h' or 'snaphu --info' to decide
% whether SNAPHU is installed: SNAPHU intentionally returns a non-zero status
% for these informational modes in current releases. Instead, first resolve
% the executable from PATH and then parse the version string from --info.
[path_status, snaphu_path] = system('command -v snaphu 2>/dev/null');
snaphu_path = strtrim(snaphu_path);

if path_status ~= 0 || isempty(snaphu_path)
    error('StaMPS_HPC:SnaphuNotFound', ...
          ['Error: "snaphu" command not found in PATH. ', ...
           'Install SNAPHU or activate the environment that provides it.']);
end

[~, cmdout] = system('snaphu --info 2>&1');
version_tokens = regexpi(cmdout, ...
    'snaphu\s+v(\d+)\.(\d+)(?:\.(\d+))?', 'tokens', 'once');

if isempty(version_tokens)
    error('StaMPS_HPC:SnaphuVersionUnknown', ...
          'Could not parse the SNAPHU version from: %s', strtrim(cmdout));
end

snaphu_major = str2double(version_tokens{1});
snaphu_minor = str2double(version_tokens{2});
if numel(version_tokens) >= 3 && ~isempty(version_tokens{3})
    snaphu_patch = str2double(version_tokens{3});
else
    snaphu_patch = 0;
end

fprintf('Detected SNAPHU v%d.%d.%d: %s\n', ...
        snaphu_major, snaphu_minor, snaphu_patch, snaphu_path);

% '-S' (single-tile re-optimization after multi-tile initialization) is a
% SNAPHU 2.x feature. Because this implementation intentionally retains -S,
% SNAPHU < 2.0.0 is not supported.
is_pre_v2 = snaphu_major < 2;
if is_pre_v2
    error('StaMPS_HPC:OldSnaphu', ...
          ['Detected SNAPHU v%d.%d.%d. This StaMPS-HPC configuration ', ...
           'retains the "-S" option and therefore requires SNAPHU >= 2.0.0.'], ...
          snaphu_major, snaphu_minor, snaphu_patch);
end

% Production recommendation:
%   2.0.0-2.0.5 : compatible with -S, but not recommended for production.
%   2.0.6       : allowed; includes the optimizer-loop fix, but 2.0.7+ is
%                 preferred for the tiled workflow.
%   2.0.7+      : recommended.
is_v2_0_pre_6 = (snaphu_major == 2 && snaphu_minor == 0 && snaphu_patch < 6);
is_v2_0_6 = (snaphu_major == 2 && snaphu_minor == 0 && snaphu_patch == 6);

if is_v2_0_pre_6
    warning('StaMPS_HPC:SnaphuLegacyProductionVersion', ...
        ['Detected SNAPHU v%d.%d.%d. The version is compatible with the ', ...
         'current -S workflow, but versions <=2.0.5 are not recommended ', ...
         'for production because later releases include important optimizer ', ...
         'loop fixes. SNAPHU v2.0.7+ is recommended.'], ...
        snaphu_major, snaphu_minor, snaphu_patch);
elseif is_v2_0_6
    warning('StaMPS_HPC:SnaphuUpgradeRecommended', ...
        ['Detected SNAPHU v2.0.6. This version is allowed, but SNAPHU ', ...
         'v2.0.7+ is recommended for the current tiled production workflow.']);
end

% --- Internal Parameters (preserved from StaMPS formulation) ---
costscale = 100;      % Scaling factor for integer cost conversion
nshortcycle = 200;    % Weighting factor for phase jumps
maxshort = 32000;     % Cost ceiling (int16 range safety)

% --- SNAPHU Tile Parameters ---
% Tile mode is enabled by NTILEROW/NTILECOL below, NOT by '-S'.
% NPROC=1 is intentional because parfor already distributes IFGs across
% MATLAB workers. This avoids nested parallel process oversubscription.
tiles_r = 3;
tiles_c = 3;
row_overlap = 100;
col_overlap = 100;
snaphu_nproc = 1;

fprintf('Unwrapping in space (HPC Optimized Mode)...\n')

% =========================================================================
% LOAD DATA
% =========================================================================
uw = load('uw_grid', 'ph', 'nzix', 'pix_size', 'n_ps', 'n_ifg');
ui = load('uw_interp');
ut = load('uw_space_time', 'dph_space_uw', 'dph_noise', 'spread', 'predef_ix');

if isfield(options, 'subset_ifg_index') && ~isempty(options.subset_ifg_index)
    subset_ifg_index = options.subset_ifg_index;
else
    subset_ifg_index = 1:size(uw.ph, 2);
end
subset_ifg_index = subset_ifg_index(:).';

if isempty(subset_ifg_index) || ...
        any(~isfinite(subset_ifg_index)) || ...
        any(subset_ifg_index ~= round(subset_ifg_index)) || ...
        any(subset_ifg_index < 1) || ...
        any(subset_ifg_index > uw.n_ifg) || ...
        numel(unique(subset_ifg_index)) ~= numel(subset_ifg_index)
    error('StaMPS_HPC:InvalidSubset', ...
          'subset_ifg_index must contain unique integer IFG indices in [1, n_ifg].');
end

if isfield(ut, 'predef_ix') && ~isempty(ut.predef_ix)
    predef_ix = ut.predef_ix;
    predef_flag = 'y';
else
    predef_flag = 'n';
    predef_ix = [];
end

[nrow, ncol] = size(uw.nzix);
disp(['Image Size: nrow = ', num2str(nrow), ', ncol = ', num2str(ncol)])

% =========================================================================
% PREPARE COST GRIDS
% =========================================================================
[y, x] = find(uw.nzix);
Z = ui.Z;
colix = ui.colix;
rowix = ui.rowix;

grid_edges = [colix(abs(colix) > 0); rowix(abs(rowix) > 0)];
n_edges = hist(abs(grid_edges), 1:ui.n_edge)';

% --- Calculate Noise/Variance ---
if strcmpi(unwrap_method, '2D')
    edge_length = sqrt(diff(x(ui.edgs(:, 2:3)), [], 2).^2 + ...
                       diff(y(ui.edgs(:, 2:3)), [], 2).^2);
    if uw.pix_size == 0
        pix_size = 5;
    else
        pix_size = uw.pix_size;
    end

    if isempty(variance)
        sigsq_noise = zeros(size(edge_length));
    else
        sigsq_noise = variance(ui.edgs(:, 2)) + variance(ui.edgs(:, 3));
    end

    sigsq_aps = (2*pi)^2;
    aps_range = 20000;
    sigsq_noise = sigsq_noise + ...
        sigsq_aps * (1 - exp(-edge_length * pix_size * 3 / aps_range));
    sigsq_noise = sigsq_noise / 10;
    dph_smooth = ut.dph_space_uw;
else
    sigsq_noise = (std(ut.dph_noise, 0, 2) / 2 / pi).^2;
    dph_smooth = ut.dph_space_uw - ut.dph_noise;
end
ut = rmfield(ut, 'dph_noise'); % Free memory

% === Optimization: Vectorized Bad Node Removal ===
is_bad_node = isnan(sigsq_noise);

valid_mask_r = ~isnan(rowix) & rowix ~= 0;
ids_r = abs(rowix(valid_mask_r));
bad_mask_r = is_bad_node(ids_r);
rowix_temp = rowix(valid_mask_r);
rowix_temp(bad_mask_r) = NaN;
rowix(valid_mask_r) = rowix_temp;

valid_mask_c = ~isnan(colix) & colix ~= 0;
ids_c = abs(colix(valid_mask_c));
bad_mask_c = is_bad_node(ids_c);
colix_temp = colix(valid_mask_c);
colix_temp(bad_mask_c) = NaN;
colix(valid_mask_c) = colix_temp;
% === End Optimization ===

sigsq = int16(round((sigsq_noise * nshortcycle^2) / costscale .* n_edges));
sigsq(sigsq < 1) = 1;

% --- Create Template Cost Matrices ---
rowcost_tmpl = zeros(nrow - 1, ncol * 4, 'int16');
colcost_tmpl = zeros(nrow, (ncol - 1) * 4, 'int16');

nzrowix = abs(rowix) > 0;
nzcolix = abs(colix) > 0;

rowcost_tmpl(:, 3:4:end) = maxshort;
colcost_tmpl(:, 3:4:end) = maxshort;

stats_ix = ~isnan(rowix);
rowcost_tmpl(:, 4:4:end) = int16(stats_ix) * (-1 - maxshort) + 1;
stats_ix = ~isnan(colix);
colcost_tmpl(:, 4:4:end) = int16(stats_ix) * (-1 - maxshort) + 1;

% --- Pre-allocate Sliced Output ---
num_total = length(subset_ifg_index);
ph_uw_slice = zeros(uw.n_ps, num_total, 'single');
msd_slice = zeros(num_total, 1);
failed_slice = false(num_total, 1);

% --- Pointers for Parfor Efficiency (broadcast/sliced variables) ---
uw_ph_sliced = uw.ph(:, subset_ifg_index);
spread_sliced = ut.spread(:, subset_ifg_index);
dph_smooth_sliced = dph_smooth(:, subset_ifg_index);
dph_space_uw_sliced = ut.dph_space_uw(:, subset_ifg_index);
uw_nzix_ptr = uw.nzix;

% =========================================================================
% PARALLEL PROCESSING LOOP
% =========================================================================
pool = gcp('nocreate');
pool_created_here = isempty(pool);
if pool_created_here
    desired_workers = feature('numCores');
    pool = parpool('local', desired_workers);
end
fprintf('Parallel processing on %d MATLAB workers...\n', pool.NumWorkers);

% --- HPC Parallel Progress Monitor Setup ---
q = parallel.pool.DataQueue;
monitor_handle = hpc_log_progress(num_total, 10, 'SNAPHU');
afterEach(q, monitor_handle);

fprintf('Starting parallel unwrapping for %d interferograms...\n', num_total);

parfor idx = 1:num_total
    i1 = subset_ifg_index(idx);
    job_ok = false;
    job_dir = sprintf('snaphu_job_%d', i1);

    try
        % 1. Create a clean isolated workspace per IFG.
        % Remove stale content from an interrupted earlier run so old output
        % cannot be mistaken for the result of the current SNAPHU execution.
        if exist(job_dir, 'dir')
            [rm_ok, rm_msg] = rmdir(job_dir, 's');
            if ~rm_ok
                error('StaMPS_HPC:JobDirCleanupFailed', ...
                      'Cannot remove stale job directory %s: %s', job_dir, rm_msg);
            end
        end

        [mk_ok, mk_msg] = mkdir(job_dir);
        if ~mk_ok
            error('StaMPS_HPC:JobDirCreateFailed', ...
                  'Cannot create %s: %s', job_dir, mk_msg);
        end

        f_conf = fullfile(job_dir, 'snaphu.conf');
        f_in = fullfile(job_dir, 'snaphu.in');
        f_out = fullfile(job_dir, 'snaphu.out');
        f_cost = fullfile(job_dir, 'snaphu.costinfile');
        f_log = fullfile(job_dir, 'snaphu.log');

        % 2. Build cost files specific to current IFG.
        rowcost = rowcost_tmpl;
        colcost = colcost_tmpl;

        spread_val = full(spread_sliced(:, idx));
        spread_val = int16(round((abs(spread_val) * nshortcycle^2) / 6 / costscale .* ...
            repmat(n_edges, 1, size(spread_val, 2))));
        sigsqtot = sigsq + spread_val;

        if predef_flag == 'y'
            sigsqtot(predef_ix(:, i1)) = 1;
        end

        rowstdgrid = ones(size(rowix), 'int16');
        rowstdgrid(nzrowix) = sigsqtot(abs(rowix(nzrowix)));
        rowcost(:, 2:4:end) = rowstdgrid;

        colstdgrid = ones(size(colix), 'int16');
        colstdgrid(nzcolix) = sigsqtot(abs(colix(nzcolix)));
        colcost(:, 2:4:end) = colstdgrid;

        % Calculate offset cycle based on smoothed phase.
        offset_cycle = (angle(exp(1i * dph_space_uw_sliced(:, idx))) - ...
                        dph_smooth_sliced(:, idx)) / 2 / pi;

        offgrid = zeros(size(rowix), 'int16');
        offgrid(nzrowix) = round(offset_cycle(abs(rowix(nzrowix))) .* ...
            sign(rowix(nzrowix)) * nshortcycle);
        rowcost(:, 1:4:end) = -offgrid;

        offgrid = zeros(size(colix), 'int16');
        offgrid(nzcolix) = round(offset_cycle(abs(colix(nzcolix))) .* ...
            sign(colix(nzcolix)) * nshortcycle);
        colcost(:, 1:4:end) = offgrid;

        % Write SNAPHU cost file using the platform's native byte order.
        fid_cost = fopen(f_cost, 'w', 'n');
        if fid_cost < 0
            error('StaMPS_HPC:CostFileOpenFailed', ...
                  'Cannot open cost file for IFG %d: %s', i1, f_cost);
        end
        fwrite(fid_cost, rowcost', 'int16');
        fwrite(fid_cost, colcost', 'int16');
        fclose(fid_cost);

        % 3. Data preparation (parallel-safe inline equivalent of writecpx).
        col_data = uw_ph_sliced(:, idx);
        ifgw = reshape(col_data(Z), nrow, ncol);

        % Transpose without conjugation so MATLAB column-major writing gives
        % the row-major sample order expected by SNAPHU.
        ifgw_t = ifgw.';

        val_to_write = zeros(2, numel(ifgw_t), 'single');
        val_to_write(1, :) = real(ifgw_t(:));
        val_to_write(2, :) = imag(ifgw_t(:));

        fid_in_handle = fopen(f_in, 'w', 'n');
        if fid_in_handle < 0
            error('StaMPS_HPC:InputFileOpenFailed', ...
                  'Cannot open input file for IFG %d: %s', i1, f_in);
        end
        fwrite(fid_in_handle, val_to_write, 'single');
        fclose(fid_in_handle);

        % 4. Generate SNAPHU configuration.
        fid_conf = fopen(f_conf, 'w');
        if fid_conf < 0
            error('StaMPS_HPC:ConfigFileOpenFailed', ...
                  'Cannot open configuration file for IFG %d: %s', i1, f_conf);
        end

        fprintf(fid_conf, 'INFILE  %s\n', f_in);
        fprintf(fid_conf, 'OUTFILE %s\n', f_out);
        fprintf(fid_conf, 'COSTINFILE %s\n', f_cost);
        fprintf(fid_conf, 'INFILEFORMAT  COMPLEX_DATA\n');
        fprintf(fid_conf, 'OUTFILEFORMAT FLOAT_DATA\n');
        fprintf(fid_conf, 'STATCOSTMODE  DEFO\n');

        % Historical StaMPS-compatible deformation parameter.
        % DEFOMAX_CYCLE controls the assumed maximum deformation-phase jump
        % when SNAPHU constructs deformation statistical costs internally.
        % Here StaMPS supplies precomputed row/column costs through COSTINFILE,
        % so this parameter normally does not rebuild those costs; it is kept
        % explicitly for compatibility and reproducibility with the inherited
        % processing configuration.
        fprintf(fid_conf, 'DEFOMAX_CYCLE %.1f\n', 0.0);

        % Tile mode configuration. These parameters, rather than '-S', enable
        % the 3x3 tiled initialization used by the current HPC workflow.
        fprintf(fid_conf, 'NTILEROW  %d\n', tiles_r);
        fprintf(fid_conf, 'NTILECOL  %d\n', tiles_c);
        fprintf(fid_conf, 'ROWOVRLP  %d\n', row_overlap);
        fprintf(fid_conf, 'COLOVRLP  %d\n', col_overlap);
        fprintf(fid_conf, 'NPROC     %d\n', snaphu_nproc);
        fclose(fid_conf);

        % 5. Execute SNAPHU.
        % IMPORTANT: '-S' is intentionally retained in v1.1.
        % It does NOT switch on tile mode. With NTILEROW/NTILECOL > 1, SNAPHU
        % first obtains the multi-tile solution and then '-S' performs an
        % additional whole-image single-tile re-optimization initialized from
        % that tiled solution. This may improve/refine the globally optimized
        % solution but also adds whole-image optimization cost. A future
        % StaMPS-HPC version may expose this behavior through an options field.
        % '-d' is retained because the externally supplied COSTINFILE does not
        % itself encode the statistical-cost mode.
        cmdstr = sprintf('snaphu -S -d -f %s %d > %s 2>&1', ...
                         f_conf, ncol, f_log);
        [snaphu_status, ~] = system(cmdstr);

        if snaphu_status ~= 0
            error('StaMPS_HPC:SnaphuExecutionFailed', ...
                  'SNAPHU failed for IFG %d with status %d. Check %s', ...
                  i1, snaphu_status, f_log);
        end

        % 6. Read and validate output.
        fid_out_handle = fopen(f_out, 'r', 'n');
        if fid_out_handle < 0
            error('StaMPS_HPC:OutputFileOpenFailed', ...
                  'Cannot read SNAPHU output for IFG %d: %s', i1, f_out);
        end

        ifguw_raw = fread(fid_out_handle, [ncol, inf], 'float32');
        fclose(fid_out_handle);

        if size(ifguw_raw, 1) ~= ncol || size(ifguw_raw, 2) ~= nrow
            error('StaMPS_HPC:InvalidSnaphuOutputSize', ...
                  ['Invalid SNAPHU output size for IFG %d: got %dx%d, ', ...
                   'expected %dx%d. Check %s'], ...
                  i1, size(ifguw_raw, 2), size(ifguw_raw, 1), ...
                  nrow, ncol, f_log);
        end

        ifguw_raw = ifguw_raw';

        % MSD calculation (preserved from original StaMPS logic).
        ifg_diff1 = ifguw_raw(1:end-1, :) - ifguw_raw(2:end, :);
        ifg_diff1 = ifg_diff1(ifg_diff1 ~= 0);
        ifg_diff2 = ifguw_raw(:, 1:end-1) - ifguw_raw(:, 2:end);
        ifg_diff2 = ifg_diff2(ifg_diff2 ~= 0);
        n_diff = length(ifg_diff1) + length(ifg_diff2);

        if n_diff > 0
            msd_slice(idx) = ...
                (sum(ifg_diff1.^2) + sum(ifg_diff2.^2)) / n_diff;
        else
            msd_slice(idx) = NaN;
        end

        % Map unwrapped grid back to PS vector.
        ph_uw_slice(:, idx) = ifguw_raw(uw_nzix_ptr);
        job_ok = true;

    catch ME
        fprintf(2, 'IFG %d failed in uw_stat_costs: %s\n', i1, ME.message);

        % Keep a compact MATLAB-side diagnostic next to the SNAPHU log.
        try
            if ~exist(job_dir, 'dir')
                mkdir(job_dir);
            end
            fid_err = fopen(fullfile(job_dir, 'matlab_error.log'), 'w');
            if fid_err >= 0
                fprintf(fid_err, 'Identifier: %s\n', ME.identifier);
                fprintf(fid_err, 'Message: %s\n', ME.message);
                fclose(fid_err);
            end
        catch
            % Diagnostic logging must never mask the original failure.
        end
    end

    failed_slice(idx) = ~job_ok;

    if ~job_ok
        ph_uw_slice(:, idx) = 0;
        msd_slice(idx) = NaN;
        % Keep failed snaphu_job_* directory and logs for diagnosis.
    else
        try
            rmdir(job_dir, 's');
        catch
            % Cleanup failure must not invalidate a scientifically valid result.
        end
    end

    send(q, 1);
end

% Close only a pool created by this function.
if pool_created_here && ~isempty(pool)
    delete(pool);
end

% Never silently save a mixture of valid and failed interferograms.
if any(failed_slice)
    failed_ifg = subset_ifg_index(failed_slice);
    error('StaMPS_HPC:SnaphuFailed', ...
          ['SNAPHU failed or produced invalid output for IFG(s): %s. ', ...
           'uw_phaseuw.mat was NOT saved. Failed snaphu_job_* directories ', ...
           'were retained for diagnosis.'], mat2str(failed_ifg));
end

% =========================================================================
% REASSEMBLE OUTPUT
% =========================================================================
ph_uw = zeros(uw.n_ps, uw.n_ifg, 'single');
msd = zeros(uw.n_ifg, 1);

for idx = 1:num_total
    ifg_id = subset_ifg_index(idx);
    ph_uw(:, ifg_id) = ph_uw_slice(:, idx);
    msd(ifg_id) = msd_slice(idx);
end

save('uw_phaseuw', 'ph_uw', 'msd')
fprintf('Unwrapping finished.\n');
end
