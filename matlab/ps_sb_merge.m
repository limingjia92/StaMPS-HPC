function [] = ps_sb_merge(ps_wd, sb_wd, merged_wd)
%PS_SB_MERGE Merge PS and SB pixels into a combined MTI dataset.
%
%   Usage:
%       ps_sb_merge(ps_wd, sb_wd)
%       ps_sb_merge(ps_wd, sb_wd, merged_wd)
%   
%   Inputs:
%       ps_wd     - Required. Path to the PS processed directory.
%       sb_wd     - Required. Path to the SBAS processed directory.
%       merged_wd - Optional. Path to save merged results. Defaults to 
%                   './MERGED' in the current working directory.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          August 16, 2026
%   Version:       1.1.0
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Memory Efficiency: Replaced dense 2D grid mapping with `ismember` 
%      for spatial intersection, drastically reducing RAM usage.
%   2. Compute Efficiency: Utilized sparse design matrices (G) and vectorized
%      operations to accelerate phase projection and SNR-weighted merging.
%   3. Data Integrity: Added precision safety for sparse operations and 
%      incorporated merging logic for raw phase datasets (ph2.mat).
%   4. Path Robustness: Replaced OS shell commands with native MATLAB I/O 
%      and implemented context-safe parameter retrieval.
%
%   Correctness & Safety Fixes (v1.1.0):
%   5. MTI pm2 Integrity: Save only the merged ph_patch field, matching the
%      original StaMPS ps_sb_merge behavior and preventing SB-only K_ps,
%      C_ps, ph_res, and coh_ps arrays from leaking into the MERGED dataset.
%   6. SB Parameter Enforcement: Require and copy the SBAS parms.mat only;
%      abort if small_baseline_flag is not 'y'.
%   7. Input Validation: Add hard checks for acquisition dates, master date,
%      image/IFG counts, IFG design indices, and key matrix dimensions.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, May 2007
%   ======================================================================

logit;

% 1. Input Parsing and Directory Setup
if nargin < 2
    error('ps_sb_merge requires at least ps_wd and sb_wd as inputs.');
end
if nargin < 3 || isempty(merged_wd)
    merged_wd = fullfile(pwd, 'MERGED');
end

if ~exist(ps_wd, 'dir')
    error('StaMPS-HPC: PS directory does not exist: %s', ps_wd);
end
if ~exist(sb_wd, 'dir')
    error('StaMPS-HPC: SBAS directory does not exist: %s', sb_wd);
end

required_input_files = {'ps2.mat', 'pm2.mat', 'rc2.mat', 'bp2.mat'};
for k = 1:numel(required_input_files)
    if ~exist(fullfile(ps_wd, required_input_files{k}), 'file')
        error('StaMPS-HPC: Missing PS input file: %s', ...
              fullfile(ps_wd, required_input_files{k}));
    end
    if ~exist(fullfile(sb_wd, required_input_files{k}), 'file')
        error('StaMPS-HPC: Missing SBAS input file: %s', ...
              fullfile(sb_wd, required_input_files{k}));
    end
end

% MERGED is an SB-network dataset. Its parameters must therefore come from
% the SBAS branch; falling back to PS parms.mat would make Step 6 follow the
% wrong processing branch.
sb_parms_file = fullfile(sb_wd, 'parms.mat');
if ~exist(sb_parms_file, 'file')
    error('StaMPS-HPC: SBAS parms.mat is required for MTI merge: %s', sb_parms_file);
end

sb_parms = load(sb_parms_file);
if ~isfield(sb_parms, 'small_baseline_flag') || ...
        ~strcmpi(strtrim(char(sb_parms.small_baseline_flag)), 'y')
    error(['StaMPS-HPC: SBAS parms.mat must contain ', ...
           'small_baseline_flag = ''y'' for MTI processing.']);
end

if ~exist(merged_wd, 'dir')
    mkdir(merged_wd);
end

[copy_ok, copy_msg] = copyfile(sb_parms_file, fullfile(merged_wd, 'parms.mat'));
if ~copy_ok
    error('StaMPS-HPC: Failed to copy SBAS parms.mat: %s', copy_msg);
end
clear sb_parms;

psver = 2;
save(fullfile(merged_wd, 'psver.mat'), 'psver', '-v7.3');

% 2. Load Core Data and Validate PS/SB Network Consistency
fprintf('Loading PS and SB core data...\n');
psps = load(fullfile(ps_wd, 'ps2.mat'));
pssb = load(fullfile(sb_wd, 'ps2.mat'));

% Required core fields. Fail early rather than allowing a mismatched PS/SB
% stack to produce a numerically valid but physically incorrect MTI result.
required_ps_fields = {'n_ps', 'n_image', 'master_ix', 'master_day', ...
                      'day', 'ij', 'lonlat'};
required_sb_fields = {'n_ps', 'n_image', 'n_ifg', 'master_day', ...
                      'day', 'ij', 'lonlat', 'ifgday_ix'};

for k = 1:numel(required_ps_fields)
    if ~isfield(psps, required_ps_fields{k})
        error('StaMPS-HPC: PS ps2.mat missing required field: %s', required_ps_fields{k});
    end
end
for k = 1:numel(required_sb_fields)
    if ~isfield(pssb, required_sb_fields{k})
        error('StaMPS-HPC: SBAS ps2.mat missing required field: %s', required_sb_fields{k});
    end
end

% Point-array dimensions.
if size(psps.ij, 1) ~= psps.n_ps || size(psps.ij, 2) < 3
    error('StaMPS-HPC: PS ij dimensions are inconsistent with n_ps.');
end
if size(pssb.ij, 1) ~= pssb.n_ps || size(pssb.ij, 2) < 3
    error('StaMPS-HPC: SBAS ij dimensions are inconsistent with n_ps.');
end
if size(psps.lonlat, 1) ~= psps.n_ps || size(psps.lonlat, 2) < 2
    error('StaMPS-HPC: PS lonlat dimensions are inconsistent with n_ps.');
end
if size(pssb.lonlat, 1) ~= pssb.n_ps || size(pssb.lonlat, 2) < 2
    error('StaMPS-HPC: SBAS lonlat dimensions are inconsistent with n_ps.');
end

% Acquisition dates and master must be identical between the two branches.
if numel(psps.day) ~= psps.n_image || numel(pssb.day) ~= pssb.n_image
    error('StaMPS-HPC: day vector length is inconsistent with n_image.');
end
if psps.n_image ~= pssb.n_image
    error('StaMPS-HPC: PS/SBAS n_image mismatch (%d vs %d).', ...
          psps.n_image, pssb.n_image);
end
if ~isequal(psps.day(:), pssb.day(:))
    error('StaMPS-HPC: PS and SBAS acquisition-date vectors are not identical.');
end

if ~isscalar(psps.master_ix) || psps.master_ix ~= round(psps.master_ix) || ...
        psps.master_ix < 1 || psps.master_ix > psps.n_image
    error('StaMPS-HPC: PS master_ix is invalid.');
end
ps_master_day = psps.day(psps.master_ix);
if ~isequal(psps.master_day, ps_master_day)
    error('StaMPS-HPC: PS master_day does not match day(master_ix).');
end
if ~isequal(pssb.master_day, ps_master_day)
    error('StaMPS-HPC: PS and SBAS master_day values are not identical.');
end
if isfield(pssb, 'master_ix')
    if ~isscalar(pssb.master_ix) || pssb.master_ix ~= psps.master_ix
        error('StaMPS-HPC: PS and SBAS master_ix values are inconsistent.');
    end
end

% SB interferogram graph dimensions and index validity.
if ~isscalar(pssb.n_ifg) || pssb.n_ifg < 1 || pssb.n_ifg ~= round(pssb.n_ifg)
    error('StaMPS-HPC: SBAS n_ifg is invalid.');
end
if size(pssb.ifgday_ix, 1) ~= pssb.n_ifg || size(pssb.ifgday_ix, 2) ~= 2
    error('StaMPS-HPC: SBAS ifgday_ix must be n_ifg-by-2.');
end
if any(~isfinite(pssb.ifgday_ix(:))) || ...
        any(pssb.ifgday_ix(:) ~= round(pssb.ifgday_ix(:))) || ...
        any(pssb.ifgday_ix(:) < 1) || any(pssb.ifgday_ix(:) > pssb.n_image)
    error('StaMPS-HPC: SBAS ifgday_ix contains invalid acquisition indices.');
end
if any(pssb.ifgday_ix(:,1) == pssb.ifgday_ix(:,2))
    error('StaMPS-HPC: SBAS ifgday_ix contains a self-interferogram.');
end

fprintf(['Input validation passed: %d acquisitions, %d SB interferograms, ', ...
         '%d PS pixels, %d SB pixels.\n'], ...
        psps.n_image, pssb.n_ifg, psps.n_ps, pssb.n_ps);

% Map spatial coordinates
[Lia, Locb] = ismember(psps.ij(:, 2:3), pssb.ij(:, 2:3), 'rows');

psu_ix  = find(~Lia);        % PS Unique
psnu_ix = find(Lia);         % PS Non-Unique (Overlapping)
sbnu_ix = Locb(Lia);         % SB Non-Unique (Overlapping)

% Resolve DORIS duplicates
[dups, I] = setdiff(psps.lonlat(psu_ix,:), pssb.lonlat, 'rows');
fprintf('%d pixels with duplicate lon/lat dropped\n', length(psu_ix) - length(dups));
psu_ix = psu_ix(I);

% 3. Construct Sparse Design Matrix (G)
i_idx = [1:pssb.n_ifg, 1:pssb.n_ifg]';
j_idx = [pssb.ifgday_ix(:,1); pssb.ifgday_ix(:,2)];
v_idx = [-ones(pssb.n_ifg, 1); ones(pssb.n_ifg, 1)];
G = sparse(i_idx, j_idx, v_idx, pssb.n_ifg, pssb.n_image);

% 4. Process pm2.mat (Phase & Coherence)
fprintf('Processing pm2.mat...\n');
pmps = load(fullfile(ps_wd, 'pm2.mat'));
pmsb = load(fullfile(sb_wd, 'pm2.mat'));

required_pmps_fields = {'ph_res', 'C_ps', 'ph_patch'};
required_pmsb_fields = {'coh_ps', 'ph_patch'};
for k = 1:numel(required_pmps_fields)
    if ~isfield(pmps, required_pmps_fields{k})
        error('StaMPS-HPC: PS pm2.mat missing required field: %s', required_pmps_fields{k});
    end
end
for k = 1:numel(required_pmsb_fields)
    if ~isfield(pmsb, required_pmsb_fields{k})
        error('StaMPS-HPC: SBAS pm2.mat missing required field: %s', required_pmsb_fields{k});
    end
end

if size(pmps.ph_res, 1) ~= psps.n_ps || size(pmps.ph_res, 2) ~= psps.n_image - 1
    error('StaMPS-HPC: PS pm2.ph_res must be n_ps-by-(n_image-1).');
end
if size(pmps.C_ps, 1) ~= psps.n_ps || size(pmps.C_ps, 2) ~= 1
    error('StaMPS-HPC: PS pm2.C_ps must be n_ps-by-1.');
end
if size(pmps.ph_patch, 1) ~= psps.n_ps || ...
        size(pmps.ph_patch, 2) ~= psps.n_image - 1
    error('StaMPS-HPC: PS pm2.ph_patch must be n_ps-by-(n_image-1).');
end
if numel(pmsb.coh_ps) ~= pssb.n_ps
    error('StaMPS-HPC: SBAS pm2.coh_ps length is inconsistent with n_ps.');
end
if size(pmsb.ph_patch, 1) ~= pssb.n_ps || ...
        size(pmsb.ph_patch, 2) ~= pssb.n_ifg
    error('StaMPS-HPC: SBAS pm2.ph_patch must be n_ps-by-n_ifg.');
end

% Project PS residuals to SB network
pmps.ph_res_all = [pmps.ph_res(:, 1:psps.master_ix-1), pmps.C_ps, pmps.ph_res(:, psps.master_ix:end)];
pmps.ph_res3 = exp(1j * (G * pmps.ph_res_all')'); 
pmps.coh_ps2 = abs(sum(pmps.ph_res3, 2)) / pssb.n_ifg; 

% Calculate SNR weights (used for all subsequent phase merges)
psnu_coh = pmps.coh_ps2(psnu_ix);
psnu_snr = 1 ./ (1 ./ psnu_coh.^2 - 1);
sbnu_coh = pmsb.coh_ps(sbnu_ix);
sbnu_snr = 1 ./ (1 ./ sbnu_coh.^2 - 1);

% Filter unique PS pixels by SB minimum coherence threshold
ps_high_coh_ix = pmps.coh_ps2(psu_ix) > min(pmsb.coh_ps);
psu_ix = psu_ix(ps_high_coh_ix);
psu_coh = pmps.coh_ps2(psu_ix);

% Save metadata
save(fullfile(merged_wd, 'merge2.mat'), 'psu_ix', 'psu_coh', 'psnu_ix', 'sbnu_ix', 'psnu_coh', 'sbnu_coh', 'G', '-v7.3');

% 5. Update and Save ps2.mat (Coordinates)
fprintf('Processing ps2.mat...\n');
ps = pssb;
ps.ij = [psps.ij(psu_ix,:); pssb.ij];
ps.lonlat = [psps.lonlat(psu_ix,:); pssb.lonlat];
ll0 = (max(ps.lonlat) + min(ps.lonlat)) / 2;
xy = llh2local(ps.lonlat', ll0) * 1000;
xy = xy';

% Context-safe parameter retrieval
return_dir = pwd; 
try
    cd(merged_wd); 
    heading = getparm('heading');
    cd(return_dir); 
catch 
    cd(return_dir); 
    heading = 0;
end

if isempty(heading), heading = 0; end
theta = (180 - heading) * pi / 180;
if theta > pi, theta = theta - 2*pi; end

rotm = [cos(theta), sin(theta); -sin(theta), cos(theta)];
xynew = (rotm * xy')';
if max(xynew(:,1)) - min(xynew(:,1)) < max(xy(:,1)) - min(xy(:,1)) && ...
   max(xynew(:,2)) - min(xynew(:,2)) < max(xy(:,2)) - min(xy(:,2))
    xy = xynew; 
    disp(['Rotating by ', num2str(theta * 180 / pi), ' degrees']);
end

xy = single(xy);
[~, sort_ix] = sortrows(xy, [2, 1]); 
xy = xy(sort_ix, :);
xy = [[1:size(xy, 1)]', xy];
xy(:, 2:3) = round(xy(:, 2:3) * 1000) / 1000; 

ps.xy = xy;
ps.lonlat = ps.lonlat(sort_ix, :);
ps.ij = ps.ij(sort_ix, :);
ps.n_ps = size(ps.xy, 1);
save(fullfile(merged_wd, 'ps2.mat'), '-struct', 'ps', '-v7.3');

% 6. Merge pm2.mat
fprintf('Merging pm2.mat...\n');
pmps.ph_patch2 = exp(1j * G * angle([pmps.ph_patch(:, 1:psps.master_ix-1), ones(psps.n_ps,1), pmps.ph_patch(:, psps.master_ix:end)])').';
% IMPORTANT: Keep MERGED pm2.mat semantically consistent with the original
% StaMPS ps_sb_merge. Do not copy SB-only K_ps/C_ps/ph_res/coh_ps fields into
% the merged dataset because their row count/order no longer matches ps.n_ps.
pm = struct();
pm.ph_patch = pmsb.ph_patch;

tmp_patch = pmsb.ph_patch(sbnu_ix, :) .* sbnu_snr + pmps.ph_patch2(psnu_ix, :) .* psnu_snr;
mask = tmp_patch ~= 0;
tmp_patch(mask) = tmp_patch(mask) ./ abs(tmp_patch(mask));
pm.ph_patch(sbnu_ix, :) = tmp_patch;

pm.ph_patch = [pmps.ph_patch2(psu_ix, :); pm.ph_patch];
pm.ph_patch = pm.ph_patch(sort_ix, :);

if size(pm.ph_patch, 1) ~= ps.n_ps || size(pm.ph_patch, 2) ~= pssb.n_ifg
    error('StaMPS-HPC: Internal error: merged pm2.ph_patch dimensions are inconsistent.');
end
save(fullfile(merged_wd, 'pm2.mat'), '-struct', 'pm', '-v7.3');
clear pm pmps pmsb tmp_patch;

% 7. Process ph2.mat (Raw Phase)
if exist(fullfile(ps_wd, 'ph2.mat'), 'file') && exist(fullfile(sb_wd, 'ph2.mat'), 'file')
    fprintf('Processing ph2.mat...\n');
    phps = load(fullfile(ps_wd, 'ph2.mat'));
    phsb = load(fullfile(sb_wd, 'ph2.mat'));

    if ~isfield(phps, 'ph') || ~isfield(phsb, 'ph')
        error('StaMPS-HPC: ph2.mat must contain field ph in both PS and SBAS branches.');
    end
    if size(phps.ph, 1) ~= psps.n_ps
        error('StaMPS-HPC: PS ph2.ph row count is inconsistent with n_ps.');
    end
    if size(phps.ph, 2) ~= pssb.n_image - 1 && size(phps.ph, 2) ~= pssb.n_image
        error('StaMPS-HPC: PS ph2.ph must have n_image-1 or n_image columns.');
    end
    if size(phsb.ph, 1) ~= pssb.n_ps || size(phsb.ph, 2) ~= pssb.n_ifg
        error('StaMPS-HPC: SBAS ph2.ph must be n_ps-by-n_ifg.');
    end
    
    % Pad PS phase with master column if missing
    if size(phps.ph, 2) == pssb.n_image - 1
        ph_ps_full = [phps.ph(:, 1:psps.master_ix-1), ones(psps.n_ps, 1), phps.ph(:, psps.master_ix:end)];
    else
        ph_ps_full = phps.ph; 
    end
    
    % Project to SBAS network
    phps.ph2 = exp(1j * G * angle(ph_ps_full)').'; 
    
    ph = phsb;
    
    % Weighted Merge
    tmp_ph = phsb.ph(sbnu_ix, :) .* sbnu_snr + phps.ph2(psnu_ix, :) .* psnu_snr;
    mask = tmp_ph ~= 0;
    tmp_ph(mask) = tmp_ph(mask) ./ abs(tmp_ph(mask));
    ph.ph(sbnu_ix, :) = tmp_ph;
    
    % Append Unique
    ph.ph = [phps.ph2(psu_ix, :); ph.ph];
    ph.ph = ph.ph(sort_ix, :);
    
    save(fullfile(merged_wd, 'ph2.mat'), '-struct', 'ph', '-v7.3');
    clear ph phps phsb tmp_ph ph_ps_full;
end

% 8. Process rc2.mat
fprintf('Processing rc2.mat...\n');
rcps = load(fullfile(ps_wd, 'rc2.mat'));
rcsb = load(fullfile(sb_wd, 'rc2.mat'));

if ~isfield(rcps, 'ph_rc') || ~isfield(rcsb, 'ph_rc')
    error('StaMPS-HPC: rc2.mat must contain field ph_rc in both PS and SBAS branches.');
end
if size(rcps.ph_rc, 1) ~= psps.n_ps || size(rcps.ph_rc, 2) ~= psps.n_image
    error('StaMPS-HPC: PS rc2.ph_rc must be n_ps-by-n_image.');
end
if size(rcsb.ph_rc, 1) ~= pssb.n_ps || size(rcsb.ph_rc, 2) ~= pssb.n_ifg
    error('StaMPS-HPC: SBAS rc2.ph_rc must be n_ps-by-n_ifg.');
end

rcps.ph_rc2 = exp(1j * G * angle(rcps.ph_rc)').'; 
rc = rcsb;

tmp_rc = rcsb.ph_rc(sbnu_ix, :) .* sbnu_snr + rcps.ph_rc2(psnu_ix, :) .* psnu_snr;
mask = tmp_rc ~= 0;
tmp_rc(mask) = tmp_rc(mask) ./ abs(tmp_rc(mask));
rc.ph_rc(sbnu_ix, :) = tmp_rc;

rc.ph_rc = [rcps.ph_rc2(psu_ix, :); rc.ph_rc];
rc.ph_rc = rc.ph_rc(sort_ix, :);
save(fullfile(merged_wd, 'rc2.mat'), '-struct', 'rc', '-v7.3');
clear rc rcps rcsb tmp_rc;

% 9. Process bp2.mat
fprintf('Processing bp2.mat...\n');
bpps = load(fullfile(ps_wd, 'bp2.mat'));
bpsb = load(fullfile(sb_wd, 'bp2.mat'));

if ~isfield(bpps, 'bperp_mat') || ~isfield(bpsb, 'bperp_mat')
    error('StaMPS-HPC: bp2.mat must contain field bperp_mat in both PS and SBAS branches.');
end
if size(bpps.bperp_mat, 1) ~= psps.n_ps || ...
        size(bpps.bperp_mat, 2) ~= psps.n_image - 1
    error('StaMPS-HPC: PS bp2.bperp_mat must be n_ps-by-(n_image-1).');
end
if size(bpsb.bperp_mat, 1) ~= pssb.n_ps || ...
        size(bpsb.bperp_mat, 2) ~= pssb.n_ifg
    error('StaMPS-HPC: SBAS bp2.bperp_mat must be n_ps-by-n_ifg.');
end

% Convert to double for sparse matrix multiplication, then revert to single
bpps_full = double([bpps.bperp_mat(:, 1:psps.master_ix-1), zeros(psps.n_ps,1), bpps.bperp_mat(:, psps.master_ix:end)]);
bpps.bperp_mat2 = single((G * bpps_full')');

bp = bpsb;
bp.bperp_mat = [bpps.bperp_mat2(psu_ix, :); bp.bperp_mat];
bp.bperp_mat = bp.bperp_mat(sort_ix, :);
save(fullfile(merged_wd, 'bp2.mat'), '-struct', 'bp', '-v7.3');
clear bp bpps bpsb bpps_full;

% 10. Dynamically Merge Simple Auxiliary Variables
aux_vars = {'hgt2', 'la2', 'inc2', 'head2'};

for v = 1:length(aux_vars)
    var_name = aux_vars{v};
    field_name = var_name(1:end-1); 
    
    sb_file = fullfile(sb_wd, [var_name, '.mat']);
    ps_file = fullfile(ps_wd, [var_name, '.mat']);
    
    if exist(sb_file, 'file') && exist(ps_file, 'file')
        fprintf('Processing %s.mat...\n', var_name);
        var_ps = load(ps_file);
        var_sb = load(sb_file);

        if ~isfield(var_ps, field_name) || ~isfield(var_sb, field_name)
            error('StaMPS-HPC: %s.mat does not contain expected field %s.', ...
                  var_name, field_name);
        end
        if size(var_ps.(field_name), 1) ~= psps.n_ps || ...
                size(var_sb.(field_name), 1) ~= pssb.n_ps
            error('StaMPS-HPC: %s dimensions are inconsistent with PS/SB n_ps.', var_name);
        end
        
        out_var = var_sb;
        out_var.(field_name) = [var_ps.(field_name)(psu_ix, :); out_var.(field_name)];
        out_var.(field_name) = out_var.(field_name)(sort_ix, :);
        
        save(fullfile(merged_wd, [var_name, '.mat']), '-struct', 'out_var', '-v7.3');
        clear var_ps var_sb out_var;
    end
end

fprintf(['MTI integrity summary: PS input=%d, SB input=%d, exact overlap=%d, ', ...
         'unique PS retained=%d, MERGED=%d.\n'], ...
        psps.n_ps, pssb.n_ps, numel(psnu_ix), numel(psu_ix), ps.n_ps);

logit(1);
fprintf('--- MTI Merge Completed Successfully (StaMPS-HPC v1.1.0) ---\n');
end
