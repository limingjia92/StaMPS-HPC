function [] = uw_grid_wrapped(ph_in, xy_in, options)
%UW_GRID_WRAPPED Resample unwrapped phase to a grid and filter (HPC Optimized).
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   =====================================================================
%   Author:        Mingjia Li
%   Date:          Febuary 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Parallel Execution: Concurrent IFG processing via 'parfor' to replace serial loops.
%   2. Vectorized Gridding: Fast phase resampling using 'accumarray' instead of explicit loops.
%   3. Filter Optimization: Internalized 'wrap_filt' with pre-calculated windows to reduce overhead.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

fprintf('Resampling phase to grid (HPC Mode)...\n')

% --- Argument Parsing ---
if nargin < 2
    error('Error: Not enough arguments. Usage: uw_grid_wrapped(ph, xy, options)');
end

if nargin < 3
    options = struct(); % Fallback if options not provided
end

% Extract parameters from options 
% 1. Grid Size (Critical Mapping: options.grid_size -> pix_size)
if isfield(options, 'grid_size')
    pix_size = options.grid_size; 
else 
    pix_size = 200; % Default
end

% 2. Prefilter Window
if isfield(options, 'prefilt_win')
    prefilt_win = options.prefilt_win; 
else 
    prefilt_win = 32; 
end

% 3. Goldstein Filter Flag
if isfield(options, 'goldfilt_flag')
    goldfilt_flag = options.goldfilt_flag; 
else 
    goldfilt_flag = 'y'; 
end

% 4. Low-pass Filter Flag
if isfield(options, 'lowfilt_flag')
    lowfilt_flag = options.lowfilt_flag; 
else 
    lowfilt_flag = 'y'; 
end

% 5. Goldstein Alpha
if isfield(options, 'gold_alpha')
    gold_alpha = options.gold_alpha; 
else 
    gold_alpha = 0.8; 
end

% 6. Predefined Phase
if isfield(options, 'ph_uw_predef')
    ph_in_predef = options.ph_uw_predef; 
else 
    ph_in_predef = []; 
end

[n_ps,n_ifg]=size(ph_in);

% --- Fix for Parfor Slicing on Empty Inputs ---
% Ensure ph_in_predef has correct dimensions for slicing even if unused.
if isempty(ph_in_predef)
    predef_flag='n';
    ph_in_predef_ptr = zeros(1, n_ifg, 'single'); 
else
    predef_flag='y';
    ph_in_predef_ptr = ph_in_predef;
end

fprintf('   Number of interferograms  : %d\n',n_ifg);
fprintf('   Number of points per ifg  : %d\n',n_ps);

if ~isreal(ph_in) && any(ph_in(:)==0)
    error('Some phase values are zero')
end

% --- Grid Setup & grid_ij Calculation ---
% Explicitly construct grid_ij for saving to ensure compatibility with 
% downstream processes (e.g., unwrapping).
if pix_size==0
    grid_x_min=1;
    grid_y_min=1;
    n_i=max(xy_in(:,3));
    n_j=max(xy_in(:,2));
    
    grid_idx_col = xy_in(:,2);
    grid_idx_row = xy_in(:,3);
    
    grid_ij = [grid_idx_row, grid_idx_col]; 
else
    grid_x_min=min(xy_in(:,2));
    grid_y_min=min(xy_in(:,3));

    % Calculate grid indices
    grid_row = ceil((xy_in(:,3)-grid_y_min+1e-3)/pix_size);
    grid_col = ceil((xy_in(:,2)-grid_x_min+1e-3)/pix_size);
    
    % Clamp to avoid index out of bounds
    max_row = max(grid_row);
    max_col = max(grid_col);
    grid_row(grid_row==max_row) = max_row-1;
    grid_col(grid_col==max_col) = max_col-1;
    
    n_i = max(grid_row);
    n_j = max(grid_col);
    
    grid_idx_row = grid_row;
    grid_idx_col = grid_col;
    
    grid_ij = [grid_row, grid_col]; 
end

% Pre-calculate Linear Indices for fast vectorized accumulation
subs_linear = sub2ind([n_i, n_j], grid_idx_row, grid_idx_col);

if min(n_i, n_j) < prefilt_win
    error(['Grid dimension (',num2str(min(n_i, n_j)),') < prefilter window (',num2str(prefilt_win),')'])
end

% --- Parallel Processing Setup ---
% Calculate valid grid mask based on geometry
grid_count = accumarray(subs_linear, 1, [n_i*n_j, 1]);
nzix = reshape(grid_count > 0, n_i, n_j);
n_ps_grid = sum(nzix(:));

% Allocate Final Outputs
ph = zeros(n_ps_grid, n_ifg, 'single');

if strcmpi(lowfilt_flag,'y')
    ph_lowpass = zeros(n_ps_grid, n_ifg, 'single');
else
    ph_lowpass = []; 
end
if predef_flag=='y'
    ph_uw_predef = zeros(n_ps_grid, n_ifg, 'single');
else
    ph_uw_predef = [];
end

% Broadcast variables to workers
ph_in_ptr = ph_in; 

fprintf('   Parallel processing on %d workers...\n', feature('numCores'));

% =========================================================================
% PARALLEL PROCESSING LOOP
% =========================================================================
parfor i1=1:n_ifg
    % 1. Prepare Data
    if isreal(ph_in_ptr)
        ph_this = exp(1i*ph_in_ptr(:,i1));
    else
        ph_this = ph_in_ptr(:,i1);
    end 
    
    % 2. Vectorized Gridding (accumarray)
    % Replaces the slow 'for i=1:n_ps' loop with optimized C++ backend
    ph_grid_vec = accumarray(subs_linear, ph_this, [n_i*n_j, 1]);
    ph_grid = reshape(ph_grid_vec, n_i, n_j);
    
    % Handle predefined model if exists
    ph_grid_uw = [];
    if predef_flag=='y'
        ph_this_uw = ph_in_predef_ptr(:,i1);
        
        valid_uw = ~isnan(ph_this_uw);
        if any(valid_uw)
            uw_sum = accumarray(subs_linear(valid_uw), ph_this_uw(valid_uw), [n_i*n_j, 1]);
            uw_count = accumarray(subs_linear(valid_uw), 1, [n_i*n_j, 1]);
            
            mask_uw = uw_count > 0;
            grid_uw_temp = zeros(n_i*n_j, 1, 'single');
            grid_uw_temp(mask_uw) = uw_sum(mask_uw) ./ uw_count(mask_uw);
            ph_grid_uw = reshape(grid_uw_temp, n_i, n_j);
        else
            ph_grid_uw = zeros(n_i, n_j, 'single');
        end
    end

    % 3. Optimized Wrap Filter
    % Call internal function to reduce overhead
    ph_this_gold = [];
    ph_this_low = [];
    
    if strcmpi(goldfilt_flag,'y') || strcmpi(lowfilt_flag,'y')
        [ph_this_gold, ph_this_low] = wrap_filt_internal(ph_grid, prefilt_win, gold_alpha, [], lowfilt_flag);
    end
    
    % 4. Store Results (Slicing)
    if strcmpi(goldfilt_flag,'y')
        ph(:,i1) = ph_this_gold(nzix);
    else
        ph(:,i1) = ph_grid(nzix);
    end
    
    if strcmpi(lowfilt_flag,'y')
        ph_lowpass(:,i1) = ph_this_low(nzix);
    end
    
    % 5. Handle Predefined Model Diff
    if predef_flag=='y'
       ph_uw_col = ph_grid_uw(nzix);
       
       ix = ~isnan(ph_uw_col);
       ph_col = ph(:,i1);
       
       ph_diff = angle(ph_col(ix) .* conj(exp(1i*ph_uw_col(ix))));
       ph_diff(abs(ph_diff)>1) = nan; 
       
       ph_uw_col(ix) = ph_uw_col(ix) + ph_diff;
       ph_uw_predef(:,i1) = ph_uw_col;
    end
end
% =========================================================================

n_ps = n_ps_grid;
fprintf('   Number of resampled points: %d\n', n_ps);

% Reconstruct coordinate info
[nz_i, nz_j] = find(nzix);
if pix_size==0
    if size(xy_in,1) ~= n_ps
         xy = [[1:n_ps]', nz_j, nz_i]; 
    else
         xy = xy_in;
    end
else
    xy = [[1:n_ps]', (nz_j-0.5)*pix_size, (nz_i-0.5)*pix_size];
end
ij = [nz_i, nz_j];

% --- Save Results ---
fprintf('Saving results uw_grid.mat ...\n');
stamps_save('uw_grid',ph,ph_in,ph_lowpass,ph_uw_predef,ph_in_predef,xy,ij,nzix,grid_x_min,grid_y_min,n_i,n_j,n_ifg,n_ps,grid_ij,pix_size)    

end

% =========================================================================
% SUBFUNCTION: Internal Optimized Wrap Filter
% =========================================================================
function [ph_out,ph_out_low]=wrap_filt_internal(ph,n_win,alpha,n_pad,low_flag)
% Optimized version of wrap_filt:
% 1. Pre-calculates window functions outside loops to avoid redundancy.
% 2. Reduces function call overhead within the HPC environment.

    if nargin<4 || isempty(n_pad)
        n_pad=round(n_win*0.25);
    end

    if nargin<5
        low_flag='n';
    end
        
    [n_i,n_j]=size(ph);
    n_inc=floor(n_win/2);
    n_win_i=ceil(n_i/n_inc)-1;
    n_win_j=ceil(n_j/n_inc)-1;

    ph_out=zeros(size(ph), 'like', ph);
    if strcmpi(low_flag,'y')
        ph_out_low=ph_out;
    else
        ph_out_low=[];
    end
    
    % --- Pre-calculation (Moved out of loops) ---
    x=[1:n_win/2];
    [X,Y]=meshgrid(x,x);
    X=X+Y;
    wind_func_base=[X,fliplr(X)];
    wind_func_base=[wind_func_base;flipud(wind_func_base)];
    
    % Smoothing kernel
    B=gausswin(7)*gausswin(7)';
    
    % Lowpass filter kernel
    L=ifftshift(gausswin(n_win+n_pad,16)*gausswin(n_win+n_pad,16)');
    
    % Handling NaNs
    ph(isnan(ph))=0;
    
    % Pre-allocate buffer
    ph_bit = zeros(n_win+n_pad, 'like', ph);
    
    % --- Main Filtering Loop ---
    for ix1=1:n_win_i
        % Optimization: Reconstruct wind_func only for edge cases
        i1=(ix1-1)*n_inc+1;
        i2=i1+n_win-1;
        
        wf=wind_func_base;
        % Vertical Edge Handling
        if i2>n_i
            i_shift=i2-n_i;
            i2=n_i;
            i1=n_i-n_win+1;
            wf=[zeros(i_shift,n_win); wf(1:n_win-i_shift,:)];
        end
        
        for ix2=1:n_win_j
            wf2=wf;
            j1=(ix2-1)*n_inc+1;
            j2=j1+n_win-1;
            
            % Horizontal Edge Handling
            if j2>n_j
               j_shift=j2-n_j;
               j2=n_j;
               j1=n_j-n_win+1;
               wf2=[zeros(n_win,j_shift), wf2(:,1:n_win-j_shift)];
            end
            
            % Extract Patch
            ph_bit(:) = 0; % Reset buffer
            ph_bit(1:n_win,1:n_win)=ph(i1:i2,j1:j2);
            
            % FFT
            ph_fft=fft2(ph_bit);
            H=abs(ph_fft);
            
            % Smoothing (filter2)
            H=ifftshift(filter2(B,fftshift(H))); 
            
            meanH=median(H(:));
            if meanH~=0
                H=H/meanH;
            end
            H=H.^alpha;
            
            % IFFT
            ph_filt=ifft2(ph_fft.*H);
            
            % Accumulate
            ph_filt_crop = ph_filt(1:n_win,1:n_win).*wf2;
            ph_out(i1:i2,j1:j2)=ph_out(i1:i2,j1:j2)+ph_filt_crop;
            
            if strcmpi(low_flag,'y')
                ph_filt_low=ifft2(ph_fft.*L);
                ph_out_low(i1:i2,j1:j2)=ph_out_low(i1:i2,j1:j2) + ph_filt_low(1:n_win,1:n_win).*wf2;
            end
        end
    end

    ph_out=abs(ph).*exp(1i*angle(ph_out)); 
    if strcmpi(low_flag,'y')
        ph_out_low=abs(ph).*exp(1i*angle(ph_out_low)); 
    end
end