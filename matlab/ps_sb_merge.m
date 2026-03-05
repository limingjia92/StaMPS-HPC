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
%   Date:          February 2026
%   Version:       1.0 
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

if ~exist(merged_wd, 'dir')
    mkdir(merged_wd);
end

% Copy auxiliary parameters safely
if exist(fullfile(sb_wd, 'parms.mat'), 'file')
    copyfile(fullfile(sb_wd, 'parms.mat'), merged_wd); 
elseif exist(fullfile(ps_wd, 'parms.mat'), 'file') 
    copyfile(fullfile(ps_wd, 'parms.mat'), merged_wd);
else
    warning('StaMPS-HPC: parms.mat not found in either ps_wd or sb_wd.');
end

psver = 2;
save(fullfile(merged_wd, 'psver.mat'), 'psver', '-v7.3');

% 2. Load Core Data and Calculate Intersections
fprintf('Loading PS and SB core data...\n');
psps = load(fullfile(ps_wd, 'ps2.mat'));
pssb = load(fullfile(sb_wd, 'ps2.mat'));

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
pm = pmsb;

tmp_patch = pmsb.ph_patch(sbnu_ix, :) .* sbnu_snr + pmps.ph_patch2(psnu_ix, :) .* psnu_snr;
mask = tmp_patch ~= 0;
tmp_patch(mask) = tmp_patch(mask) ./ abs(tmp_patch(mask));
pm.ph_patch(sbnu_ix, :) = tmp_patch;

pm.ph_patch = [pmps.ph_patch2(psu_ix, :); pm.ph_patch];
pm.ph_patch = pm.ph_patch(sort_ix, :);
save(fullfile(merged_wd, 'pm2.mat'), '-struct', 'pm', '-v7.3');
clear pm pmps pmsb tmp_patch;

% 7. Process ph2.mat (Raw Phase)
if exist(fullfile(ps_wd, 'ph2.mat'), 'file') && exist(fullfile(sb_wd, 'ph2.mat'), 'file')
    fprintf('Processing ph2.mat...\n');
    phps = load(fullfile(ps_wd, 'ph2.mat'));
    phsb = load(fullfile(sb_wd, 'ph2.mat'));
    
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
        
        out_var = var_sb;
        out_var.(field_name) = [var_ps.(field_name)(psu_ix, :); out_var.(field_name)];
        out_var.(field_name) = out_var.(field_name)(sort_ix, :);
        
        save(fullfile(merged_wd, [var_name, '.mat']), '-struct', 'out_var', '-v7.3');
        clear var_ps var_sb out_var;
    end
end

logit(1);
fprintf('--- MTI Merge Completed Successfully ---\n');
end