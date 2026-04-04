function [] = uw_grid_wrapped(tmp_ph_file, xy_in, options, unwrap_ifg_index)
%UW_GRID_WRAPPED Resample unwrapped phase to a grid and filter (HPC Block-IO).
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   =====================================================================
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       2.0 (Block-IO & Matfile Parfor Streaming)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Two-Phase Architecture: Completely eradicated 'parfor' D-state lockups 
%      (Orthogonal Read Penalty) by decoupling disk I/O from parallel workers. 
%      Phase 1 streams Row-Blocks to RAM; Phase 2 filters in-memory.
%   2. Parfor Strict Compliance: Resolved conditional variable instantiation 
%      errors (e.g., 'predef_grid_mean') through unconditional dummy allocation.
%   3. Vectorized Gridding & Filtering: Accelerated spatial resampling via 
%      'accumarray' and internalized 'wrap_filt' to eliminate overhead.
%   4. Ghost Variable Purge: Stripped the massive 'ph_in' ghost array from 
%      the 'uw_grid.mat' save step, drastically freeing RAM and Disk I/O.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

fprintf('Resampling phase to grid (HPC Mode)...\n')

% --- Argument Parsing ---
if nargin < 2
    error('Error: Not enough arguments. Usage: uw_grid_wrapped(tmp_ph_file, xy, options, unwrap_ifg_index)');
end
if nargin < 3, options = struct(); end
if nargin < 4, error('unwrap_ifg_index is required in HPC version.'); end

% Extract parameters from options 
pix_size = 200; if isfield(options, 'grid_size'), pix_size = options.grid_size; end
prefilt_win = 32; if isfield(options, 'prefilt_win'), prefilt_win = options.prefilt_win; end
goldfilt_flag = 'y'; if isfield(options, 'goldfilt_flag'), goldfilt_flag = options.goldfilt_flag; end
lowfilt_flag = 'y'; if isfield(options, 'lowfilt_flag'), lowfilt_flag = options.lowfilt_flag; end
gold_alpha = 0.8; if isfield(options, 'gold_alpha'), gold_alpha = options.gold_alpha; end

% Predefined Phase Check
has_predef = isfield(options, 'ph_uw_predef_file') && ~isempty(options.ph_uw_predef_file);

n_ps = size(xy_in, 1);
n_ifg = length(unwrap_ifg_index);

fprintf('   Number of interferograms to process: %d\n', n_ifg);
fprintf('   Number of PS points                : %d\n', n_ps);

% --- Grid Setup & grid_ij Calculation ---
if pix_size==0
    grid_x_min=1; grid_y_min=1;
    n_i=max(xy_in(:,3)); n_j=max(xy_in(:,2));
    grid_idx_row = xy_in(:,3); grid_idx_col = xy_in(:,2);
else
    grid_x_min=min(xy_in(:,2)); grid_y_min=min(xy_in(:,3));
    grid_row = ceil((xy_in(:,3)-grid_y_min+1e-3)/pix_size);
    grid_col = ceil((xy_in(:,2)-grid_x_min+1e-3)/pix_size);
    
    max_row = max(grid_row); max_col = max(grid_col);
    grid_row(grid_row==max_row) = max_row-1;
    grid_col(grid_col==max_col) = max_col-1;
    
    n_i = max(grid_row); n_j = max(grid_col);
    grid_idx_row = grid_row; grid_idx_col = grid_col;
end
grid_ij = [grid_idx_row, grid_idx_col]; 

% Pre-calculate Linear Indices for fast vectorized accumulation
subs_linear = sub2ind([n_i, n_j], grid_idx_row, grid_idx_col);

if min(n_i, n_j) < prefilt_win
    error(['Grid dimension (',num2str(min(n_i, n_j)),') < prefilter window (',num2str(prefilt_win),')'])
end

% --- Calculate Grid Extents ---
grid_count = accumarray(subs_linear, 1, [n_i*n_j, 1]);
nzix = reshape(grid_count > 0, n_i, n_j);
n_ps_grid = sum(nzix(:));


% =========================================================================
% PHASE 1: Sequential Row-Block Gridding (Disk to RAM)
% =========================================================================
fprintf('   Phase 1/2: Sequential in-memory gridding (Row-Block Streaming)...\n');

ph_grid_sum = zeros(n_i*n_j, n_ifg, 'single');
ph_grid_sum = complex(ph_grid_sum, ph_grid_sum);

% Dummy initialization for parfor safety
predef_grid_mean = nan(n_i*n_j, n_ifg, 'single'); 
if has_predef
    predef_sum = zeros(n_i*n_j, n_ifg, 'single');
    predef_count = zeros(n_i*n_j, n_ifg, 'single');
end

m_tmp = matfile(tmp_ph_file);
row_block_size = 1000000;
n_blocks = ceil(n_ps / row_block_size);

% Calculate absolute continuous bounding box for Matfile I/O
min_ifg_idx = min(unwrap_ifg_index);
max_ifg_idx = max(unwrap_ifg_index);
relative_ifg_index = unwrap_ifg_index - min_ifg_idx + 1;

for b = 1:n_blocks
    rows = ((b-1)*row_block_size + 1) : min(b*row_block_size, n_ps);
    
    % Read contiguous block from disk, then slice irregularly in RAM
    ph_chunk_raw = m_tmp.ph_w(rows, min_ifg_idx:max_ifg_idx);
    ph_chunk = ph_chunk_raw(:, relative_ifg_index);
    clear ph_chunk_raw; 

    if isreal(ph_chunk); ph_chunk = exp(1i * ph_chunk); end
    
    subs_chunk = subs_linear(rows);
    
    for c = 1:n_ifg
        ph_grid_sum(:, c) = ph_grid_sum(:, c) + accumarray(subs_chunk, ph_chunk(:, c), [n_i*n_j, 1]);
    end
    
    if has_predef
        % Same proxy logic for predef phase
        pd_chunk_raw = m_tmp.ph_uw_predef(rows, min_ifg_idx:max_ifg_idx);
        pd_chunk = pd_chunk_raw(:, relative_ifg_index);
        clear pd_chunk_raw;

        for c = 1:n_ifg
            pd_col = pd_chunk(:, c);
            valid = ~isnan(pd_col);
            if any(valid)
                predef_sum(:, c) = predef_sum(:, c) + accumarray(subs_chunk(valid), pd_col(valid), [n_i*n_j, 1]);
                predef_count(:, c) = predef_count(:, c) + accumarray(subs_chunk(valid), 1, [n_i*n_j, 1]);
            end
        end
    end
    if mod(b,5)==0 || b==n_blocks
        fprintf('      Aggregated Block %d / %d\n', b, n_blocks);
    end
end

if has_predef
    valid_pd = predef_count > 0;
    predef_grid_mean(valid_pd) = predef_sum(valid_pd) ./ predef_count(valid_pd);
end

clear predef_sum predef_count ph_chunk pd_chunk


% =========================================================================
% PHASE 2: In-Memory Parallel Filtering (24 Cores Full Speed)
% =========================================================================
fprintf('   Phase 2/2: Parallel spatial filtering on %d workers...\n', feature('numCores'));

ph = zeros(n_ps_grid, n_ifg, 'single');
ph_lowpass_out = zeros(n_ps_grid, n_ifg, 'single');
ph_uw_predef_out = zeros(n_ps_grid, n_ifg, 'single');

% Launch parfor purely in-memory
parfor i1 = 1:n_ifg
    % 1. Extract 2D grid for this IFG (DIRECTLY USING SUM, NOT MEAN)
    ph_grid = reshape(ph_grid_sum(:, i1), n_i, n_j);
    ph_grid_uw = reshape(predef_grid_mean(:, i1), n_i, n_j); 

    % 2. Optimized Wrap Filter
    ph_this_gold = []; ph_this_low = [];
    if strcmpi(goldfilt_flag,'y') || strcmpi(lowfilt_flag,'y')
        [ph_this_gold, ph_this_low] = wrap_filt_internal(ph_grid, prefilt_win, gold_alpha, [], lowfilt_flag);
    end
    
    % 3. Store Results unconditionally to satisfy parfor
    if strcmpi(goldfilt_flag,'y')
        ph(:,i1) = ph_this_gold(nzix);
    else
        ph(:,i1) = ph_grid(nzix); 
    end
    
    if strcmpi(lowfilt_flag,'y')
        ph_lowpass_out(:,i1) = ph_this_low(nzix); 
    else
        ph_lowpass_out(:,i1) = 0; % Dummy assign
    end
    
    % 4. Handle Predefined Model Diff
    if has_predef
       ph_uw_col = ph_grid_uw(nzix);
       ix = ~isnan(ph_uw_col);
       ph_col = ph(:,i1);
       
       ph_diff = angle(ph_col(ix) .* conj(exp(1i*ph_uw_col(ix))));
       ph_diff(abs(ph_diff)>1) = nan; 
       
       ph_uw_col(ix) = ph_uw_col(ix) + ph_diff;
       ph_uw_predef_out(:,i1) = ph_uw_col;
    else
       ph_uw_predef_out(:,i1) = 0; % Dummy assign
    end
end
% =========================================================================

if strcmpi(lowfilt_flag,'y')
    ph_lowpass = ph_lowpass_out;
else
    ph_lowpass = [];
end

if has_predef
    ph_uw_predef = ph_uw_predef_out;
else
    ph_uw_predef = [];
end
clear ph_lowpass_out ph_uw_predef_out;

n_ps_orig = n_ps;
n_ps = n_ps_grid; 
fprintf('   Number of resampled points: %d\n', n_ps);

[nz_i, nz_j] = find(nzix);
if pix_size==0
    if size(xy_in,1) ~= n_ps_orig, xy = [[1:n_ps]', nz_j, nz_i]; else, xy = xy_in; end
else
    xy = [[1:n_ps]', (nz_j-0.5)*pix_size, (nz_i-0.5)*pix_size];
end
ij = [nz_i, nz_j];

% --- Save Results ---
fprintf('Saving results uw_grid.mat \n');
stamps_save('uw_grid',ph,ph_lowpass,ph_uw_predef,xy,ij,nzix,grid_x_min,grid_y_min,n_i,n_j,n_ifg,n_ps,grid_ij,pix_size);

end

% =========================================================================
% SUBFUNCTION: Internal Optimized Wrap Filter
% =========================================================================
function [ph_out,ph_out_low]=wrap_filt_internal(ph,n_win,alpha,n_pad,low_flag)
% Optimized version of wrap_filt:
% 1. Pre-calculates window functions outside loops to avoid redundancy.
% 2. Reduces function call overhead within the HPC environment.

    if nargin<4 || isempty(n_pad), n_pad=round(n_win*0.25); end
    if nargin<5, low_flag='n'; end
        
    [n_i,n_j]=size(ph);
    n_inc=floor(n_win/2);
    n_win_i=ceil(n_i/n_inc)-1;
    n_win_j=ceil(n_j/n_inc)-1;

    ph_out=zeros(size(ph), 'like', ph);
    if strcmpi(low_flag,'y'), ph_out_low=ph_out; else, ph_out_low=[]; end
    
    x=[1:n_win/2];
    [X,Y]=meshgrid(x,x);
    X=X+Y;
    wind_func_base=[X,fliplr(X)];
    wind_func_base=[wind_func_base;flipud(wind_func_base)];
    
    B=gausswin(7)*gausswin(7)';
    L=ifftshift(gausswin(n_win+n_pad,16)*gausswin(n_win+n_pad,16)');
    ph(isnan(ph))=0;
    ph_bit = zeros(n_win+n_pad, 'like', ph);
    
    for ix1=1:n_win_i
        i1=(ix1-1)*n_inc+1; i2=i1+n_win-1;
        wf=wind_func_base;
        if i2>n_i
            i_shift=i2-n_i; i2=n_i; i1=n_i-n_win+1;
            wf=[zeros(i_shift,n_win); wf(1:n_win-i_shift,:)];
        end
        for ix2=1:n_win_j
            wf2=wf; j1=(ix2-1)*n_inc+1; j2=j1+n_win-1;
            if j2>n_j
               j_shift=j2-n_j; j2=n_j; j1=n_j-n_win+1;
               wf2=[zeros(n_win,j_shift), wf2(:,1:n_win-j_shift)];
            end
            ph_bit(:) = 0; 
            ph_bit(1:n_win,1:n_win)=ph(i1:i2,j1:j2);
            ph_fft=fft2(ph_bit);
            H=abs(ph_fft);
            H=ifftshift(filter2(B,fftshift(H))); 
            meanH=median(H(:));
            if meanH~=0, H=H/meanH; end
            H=H.^alpha;
            ph_filt=ifft2(ph_fft.*H);
            ph_filt_crop = ph_filt(1:n_win,1:n_win).*wf2;
            ph_out(i1:i2,j1:j2)=ph_out(i1:i2,j1:j2)+ph_filt_crop;
            if strcmpi(low_flag,'y')
                ph_filt_low=ifft2(ph_fft.*L);
                ph_out_low(i1:i2,j1:j2)=ph_out_low(i1:i2,j1:j2) + ph_filt_low(1:n_win,1:n_win).*wf2;
            end
        end
    end
    ph_out=abs(ph).*exp(1i*angle(ph_out)); 
    if strcmpi(low_flag,'y'), ph_out_low=abs(ph).*exp(1i*angle(ph_out_low)); end
end