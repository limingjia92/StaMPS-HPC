function [] = uw_stat_costs(options)
%UW_STAT_COSTS Generate statistical costs for MAP phase unwrapping (HPC Optimized).
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   =====================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Pre-Unwrapping Filtering: Implements a weighted Goldstein-Werner filter to reduce noise
%      before unwrapping, preserving magnitude as coherence weights.
%   2. Parallel Execution: Uses parfor to process interferograms concurrently (1 IFG per Worker).
%   3. Snaphu Optimization: Leverages 'snaphu -S' (Tile Mode) with 'NPROC=1' to maximize 
%      single-threaded efficiency within the parallel framework.
%   4. Logic Optimization: Vectorized cost calculation.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, May 2007
%   ======================================================================

if nargin < 1
    options = struct(); % Fallback
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

% --- Snaphu Version Check ---
% Ensure Snaphu v2.0+ is available for the -S flag and optimized tiling.
[status, cmdout] = system('snaphu -h');
if status == 0
    error('StaMPS_HPC:SnaphuNotFound', ...
          'Error: "snaphu" command not found in system path. Please install Snaphu.');
end

if ~contains(cmdout, 'v2') && ~contains(cmdout, 'v2.0') 
    % Note: Some builds might just say "v2", others "v2.0.x"
    error('StaMPS_HPC:OldSnaphu', ...
        ['Critical Error: Detected Snaphu version < 2.0.\n' ...
         'The HPC optimization strategy relies on the "-S" flag (Tile Mode) introduced in v2.0.\n']);
end

% --- Internal Parameters ---
costscale = 100;      % Scaling factor for integer cost conversion
nshortcycle = 200;    % Weighting factor for phase jumps (smoothness constraint)
maxshort = 32000;     % Cost ceiling (int16 limit)

fprintf('Unwrapping in space (HPC Optimized Mode)...\n')

% --- Load Data ---
uw=load('uw_grid','ph','nzix','pix_size','n_ps','n_ifg');
ui=load('uw_interp');
ut=load('uw_space_time','dph_space_uw','dph_noise','spread');

if isfield(options, 'subset_ifg_index') && ~isempty(options.subset_ifg_index)
    subset_ifg_index = options.subset_ifg_index;
else
    subset_ifg_index = 1:size(uw.ph, 2);
end

if isfield(ut,'predef_ix') && ~isempty(ut.predef_ix)
    predef_ix = ut.predef_ix;
    predef_flag='y';
else
    predef_flag='n';
    predef_ix = []; 
end

[nrow,ncol]=size(uw.nzix);
disp(['Image Size: nrow = ',num2str(nrow),', ncol = ',num2str(ncol)])
fprintf('Parallel processing on %d workers...\n', feature('numCores'));

% --- Prepare Cost Grids ---
[y,x]=find(uw.nzix);
Z=ui.Z;
colix=ui.colix;
rowix=ui.rowix;

grid_edges=[colix(abs(colix)>0);rowix(abs(rowix)>0)];
n_edges=hist(abs(grid_edges),[1:ui.n_edge])';

% --- Calculate Noise/Variance ---
if strcmpi(unwrap_method,'2D')
    edge_length=sqrt(diff(x(ui.edgs(:,2:3)),[],2).^2+diff(y(ui.edgs(:,2:3)),[],2).^2);
    if uw.pix_size==0, pix_size=5; else, pix_size=uw.pix_size; end
    
    if isempty(variance)
        sigsq_noise=zeros(size(edge_length));
    else
        sigsq_noise=variance(ui.edgs(:,2))+variance(ui.edgs(:,3));
    end
    sigsq_aps=(2*pi)^2; 
    aps_range=20000; 
    sigsq_noise=sigsq_noise+sigsq_aps*(1-exp(-edge_length*pix_size*3/aps_range));
    sigsq_noise=sigsq_noise/10; 
    dph_smooth=ut.dph_space_uw;
else
    sigsq_noise=(std(ut.dph_noise,0,2)/2/pi).^2;
    dph_smooth=ut.dph_space_uw-ut.dph_noise;
end
ut=rmfield(ut,'dph_noise'); % Free memory

% === Optimization: Vectorized Bad Node Removal ===
% 1. Create Lookup Table for Bad Nodes (NaN variance)
is_bad_node = isnan(sigsq_noise); 

% 2. Batch Process Row Indices
valid_mask_r = ~isnan(rowix) & rowix ~= 0; 
ids_r = abs(rowix(valid_mask_r));
bad_mask_r = is_bad_node(ids_r);
% Apply NaN to connections involving bad nodes
rowix_temp = rowix(valid_mask_r);
rowix_temp(bad_mask_r) = NaN;
rowix(valid_mask_r) = rowix_temp;

% 3. Batch Process Column Indices
valid_mask_c = ~isnan(colix) & colix ~= 0;
ids_c = abs(colix(valid_mask_c));
bad_mask_c = is_bad_node(ids_c);
colix_temp = colix(valid_mask_c);
colix_temp(bad_mask_c) = NaN;
colix(valid_mask_c) = colix_temp;
% === End Optimization ===

sigsq=int16(round(((sigsq_noise)*nshortcycle^2)/costscale.*n_edges));
sigsq(sigsq<1)=1;

% --- Create Template Cost Matrices ---
rowcost_tmpl=zeros((nrow-1),ncol*4,'int16');
colcost_tmpl=zeros((nrow),(ncol-1)*4,'int16');

nzrowix=abs(rowix)>0;
nzcolix=abs(colix)>0;

rowcost_tmpl(:,3:4:end)= maxshort;
colcost_tmpl(:,3:4:end)= maxshort;

stats_ix=~isnan(rowix); 
rowcost_tmpl(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1; 
stats_ix=~isnan(colix); 
colcost_tmpl(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1;

% --- Filter Parameters Setup ---
% Use setparm to configure aggressive filtering for specific noisy IFGs.
uwnoisey1 = getparm('uwnoisey1', 1); 
if isempty(uwnoisey1), uwnoisey1=[]; end

uwnoisey2 = getparm('uwnoisey2', 1); 
if isempty(uwnoisey2), uwnoisey2=[]; end

% --- Pre-allocate Sliced Output ---
ph_uw_slice = zeros(uw.n_ps, length(subset_ifg_index), 'single');
msd_slice = zeros(length(subset_ifg_index), 1);

% --- Pointers for Parfor Efficiency (Broadcast variables) ---
uw_ph_ptr = uw.ph; 
uw_nzix_ptr = uw.nzix;
spread_ptr = ut.spread;

% =========================================================================
% PARALLEL PROCESSING LOOP
% =========================================================================
desired_workers = feature('numCores'); 
if isempty(gcp('nocreate'))
    parpool('local', desired_workers);
end

num_total = length(subset_ifg_index);

% --- HPC Parallel Progress Monitor Setup ---
q = parallel.pool.DataQueue;
monitor_handle = hpc_log_progress(num_total, 10, 'SNAPHU'); % print out each 10%
afterEach(q, monitor_handle); 
% -------------------------------------------

fprintf('Starting parallel unwrapping for %d interferograms...\n', num_total);

parfor idx = 1:length(subset_ifg_index)
    i1 = subset_ifg_index(idx);
    
    % 1. Create Isolated Workspace per Worker
    job_dir = sprintf('snaphu_job_%d', i1);
    if ~exist(job_dir, 'dir'), mkdir(job_dir); end
    
    f_conf = fullfile(job_dir, 'snaphu.conf');
    f_in = fullfile(job_dir, 'snaphu.in');
    f_out = fullfile(job_dir, 'snaphu.out');
    f_cost = fullfile(job_dir, 'snaphu.costinfile');
    f_log = fullfile(job_dir, 'snaphu.log');

    % 2. Build Cost Files specific to current IFG
    rowcost = rowcost_tmpl;
    colcost = colcost_tmpl;
    
    spread_val = full(spread_ptr(:,i1));
    spread_val = int16(round((abs(spread_val)*nshortcycle^2)/6/costscale.*repmat(n_edges,1,size(spread_val,2))));
    sigsqtot = sigsq + spread_val;
    
    if predef_flag=='y'
         sigsqtot(predef_ix(:,i1))=1;
    end
    
    rowstdgrid = ones(size(rowix),'int16');
    rowstdgrid(nzrowix) = sigsqtot(abs(rowix(nzrowix)));
    rowcost(:,2:4:end) = rowstdgrid;

    colstdgrid = ones(size(colix),'int16');
    colstdgrid(nzcolix) = sigsqtot(abs(colix(nzcolix)));
    colcost(:,2:4:end) = colstdgrid;
    
    % Calculate offset cycle based on smoothed phase
    offset_cycle = (angle(exp(1i*ut.dph_space_uw(:,i1))) - dph_smooth(:,i1))/2/pi;
    
    offgrid = zeros(size(rowix),'int16');
    offgrid(nzrowix) = round(offset_cycle(abs(rowix(nzrowix))).*sign(rowix(nzrowix))*nshortcycle);
    rowcost(:,1:4:end) = -offgrid;

    offgrid = zeros(size(colix),'int16');
    offgrid(nzcolix) = round(offset_cycle(abs(colix(nzcolix))).*sign(colix(nzcolix))*nshortcycle);
    colcost(:,1:4:end) = offgrid;
    
    % Write Cost File (Native Endian or Explicit LE)
    fid_cost = fopen(f_cost, 'w'); 
    fwrite(fid_cost, rowcost', 'int16');
    fwrite(fid_cost, colcost', 'int16');
    fclose(fid_cost);
    
    % 3. Filtering & Data Prep
    ifgw = reshape(uw_ph_ptr(Z,i1), nrow, ncol);
    
    % Adaptive filter settings
    alpha = 0.6; % Default alpha
    nfft = 16;   % Default FFT size
    
    if ismember(i1, uwnoisey1)
        alpha = 1.0; nfft = 32;
    elseif ismember(i1, uwnoisey2)
        alpha = 1.5; nfft = 64;
    end
    
    % Apply Weighted Goldstein-Werner Filter
    % wfrac_min set to 0.7 to handle edges
    ifgw_filt = goldstein_werner_filter_weighted(ifgw, nfft, alpha, 0.7);
    
    % CRITICAL: Transpose for Row-Major (C-style) compatibility
    ifgw_filt_t = ifgw_filt.'; 

    % Interleave Real and Imaginary parts
    val_to_write = zeros(2, numel(ifgw_filt_t), 'single');
    val_to_write(1,:) = real(ifgw_filt_t(:));
    val_to_write(2,:) = imag(ifgw_filt_t(:));
    
    % Write input file (Force Little Endian for consistency)
    fid_in_handle = fopen(f_in, 'w');
    fwrite(fid_in_handle, val_to_write, 'float32', 'ieee-le'); 
    fclose(fid_in_handle);
    
    % 4. Generate Snaphu Configuration
    fid_conf = fopen(f_conf, 'w');
    fprintf(fid_conf,'INFILE  %s\n', f_in);
    fprintf(fid_conf,'OUTFILE %s\n', f_out);
    fprintf(fid_conf,'COSTINFILE %s\n', f_cost);
    fprintf(fid_conf,'INFILEFORMAT  COMPLEX_DATA\n');
    fprintf(fid_conf,'OUTFILEFORMAT FLOAT_DATA\n');
    fprintf(fid_conf,'STATCOSTMODE  DEFO\n'); 
    fprintf(fid_conf,'DEFOMAX_CYCLE %d\n', 0.0);
    
    % --- Parallel/Tiling Parameters ---
    % Using 3x3 tiles improves performance even on single core (NPROC=1)
    tiles_r = 3; 
    tiles_c = 3;
    
    fprintf(fid_conf,'NTILEROW  %d\n', tiles_r); 
    fprintf(fid_conf,'NTILECOL  %d\n', tiles_c); 
    fprintf(fid_conf,'ROWOVRLP  %d\n', 100); % Increased overlap for safety
    fprintf(fid_conf,'COLOVRLP  %d\n', 100); 
    fprintf(fid_conf,'NPROC     1\n'); % Force single thread per worker
    fclose(fid_conf);
    
    % 5. Execute Snaphu
    % Using -S flag for tile mode optimization (Snaphu v2.0+)
    cmdstr = sprintf('snaphu -S -d -f %s %d > %s 2>&1', f_conf, ncol, f_log);
    [status, ~] = system(cmdstr);
    
    if status ~= 0
        fprintf('Error running snaphu for IFG %d. Check %s\n', i1, f_log);
        ph_uw_slice(:, idx) = 0; 
        msd_slice(idx) = NaN;
    else
        % 6. Read Output
        fid_out_handle = fopen(f_out, 'r');
        if fid_out_handle > 0
            ifguw_raw = fread(fid_out_handle, [ncol, inf], 'float32');
            fclose(fid_out_handle);
            ifguw_raw = ifguw_raw'; 
            
            % MSD Calculation
            ifg_diff1 = ifguw_raw(1:end-1,:) - ifguw_raw(2:end,:);
            ifg_diff1 = ifg_diff1(ifg_diff1~=0);
            ifg_diff2 = ifguw_raw(:,1:end-1) - ifguw_raw(:,2:end);
            ifg_diff2 = ifg_diff2(ifg_diff2~=0);
            msd_slice(idx) = (sum(ifg_diff1.^2) + sum(ifg_diff2.^2)) / (length(ifg_diff1) + length(ifg_diff2));
            
            % Map back to PS vector
            ph_uw_slice(:, idx) = ifguw_raw(uw_nzix_ptr);
        else
            fprintf('Error reading output for IFG %d\n', i1);
        end
    end
    
    % 7. Cleanup
    if status == 0
        try rmdir(job_dir, 's'); catch; end
    end

    send(q, 1); % send signal to Monitor
end
delete(gcp('nocreate'));

% --- Reassemble Output ---
ph_uw = zeros(uw.n_ps, uw.n_ifg, 'single');
msd = zeros(uw.n_ifg, 1);

for idx = 1:length(subset_ifg_index)
    ifg_id = subset_ifg_index(idx);
    ph_uw(:, ifg_id) = ph_uw_slice(:, idx);
    msd(ifg_id) = msd_slice(idx);
end

save('uw_phaseuw','ph_uw','msd')
fprintf('Unwrapping finished.\n');
end

% =========================================================================
% SUBFUNCTION: Weighted Goldstein-Werner Filter
% =========================================================================
function im_filt = goldstein_werner_filter_weighted(im, nfft, alpha, wfrac_min)
%GOLDSTEIN_WERNER_FILTER_WEIGHTED Optimized Goldstein-Werner phase filter
%
%
%   im_filt = goldstein_werner_filter_weighted(im, nfft, alpha, wfrac_min)
%
%   INPUTS:
%       im        : Complex interferogram (2D double/single)
%       nfft      : FFT window size (e.g., 32, 64)
%       alpha     : Filter exponent (0.5 = standard, >1.0 = strong)
%       wfrac_min : Minimum fraction of valid pixels in window (0.0-1.0)
%
%   OUTPUT:
%       im_filt   : Filtered complex image.
%                   NOTE: Magnitude is NOT normalized; it represents
%                   local phase quality/coherence, used as weights by Snaphu.
%
%   ALGORITHM:
%       Implements the Goldstein-Werner spectral filter with a high-overlap
%       strategy (step = nfft/8) to reduce tiling artifacts and preserve
%       coherence details.


    if nargin < 2 || isempty(nfft), nfft = 16; end % GAMMA default/usage implies 16 or 32
    if nargin < 3 || isempty(alpha), alpha = 0.6; end % GAMMA usage was 0.6
    if nargin < 4, wfrac_min = 0.7; end
    
    [rows, cols] = size(im);
    
    % ---  Goldstein Filtering  ---
    step = max(1, floor(nfft / 8)); 
    H_smooth = ones(5) / 25; 
    win_1d = triang(nfft);
    window = win_1d * win_1d';
    
    pad_r = ceil(rows/step)*step + nfft;
    pad_c = ceil(cols/step)*step + nfft;
    im_pad = zeros(pad_r, pad_c, 'like', im);
    im_pad(1:rows, 1:cols) = im;
    
    acc_im = zeros(pad_r, pad_c, 'like', im);
    acc_w = zeros(pad_r, pad_c, 'single');
    
    min_pixels = floor(wfrac_min * nfft * nfft);

    for r = 1:step:(rows)
        for c = 1:step:(cols)
            r_end = r + nfft - 1;
            c_end = c + nfft - 1;
            if r_end > pad_r || c_end > pad_c, continue; end
            
            patch = im_pad(r:r_end, c:c_end);
            if nnz(patch) < min_pixels, continue; end
            
            Z = fft2(patch);
            MagZ = abs(Z);
            S_smooth = conv2(MagZ, H_smooth, 'same');
            Z_filt = Z .* (S_smooth .^ alpha);
            patch_filt = ifft2(Z_filt);
            
            acc_im(r:r_end, c:c_end) = acc_im(r:r_end, c:c_end) + patch_filt .* window;
            acc_w(r:r_end, c:c_end) = acc_w(r:r_end, c:c_end) + window;
        end
    end
    
    mask = acc_w > 1e-6;
    acc_im(mask) = acc_im(mask) ./ acc_w(mask);
    raw_filt = acc_im(1:rows, 1:cols);
    raw_filt(isnan(raw_filt)) = 0;

    % - Magnitude Replacement  ---
    % Extract Phase
    phase_only = exp(1i * angle(raw_filt));
    
    % Calculate Coherence using a 7x7 window 
    cc_win = 7;
    kernel = ones(cc_win) / (cc_win^2);
    
    % Coherence = |mean(exp(j*phi))
    
    coh_est = abs(imfilter(phase_only, kernel, 'symmetric'));
    
    % Clamp Coherence to 0.001 - 1.0 range for Snaphu stability
    coh_est(coh_est < 0.001) = 0.001;
    coh_est(coh_est > 1.0) = 1.0;
    
    % Reconstruct Output: Phase from Filter, Magnitude from Coherence
    im_filt = coh_est .* phase_only;

end
