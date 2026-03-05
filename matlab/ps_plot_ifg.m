function [ph_lims] = ps_plot_ifg(ps, in_ph, bg_flag, col_rg, lon_rg, lat_rg)
%PS_PLOT_IFG Render PS scatter map using vectorized engine (HPC Optimized)
%   Takes directly pre-loaded PS struct and assembled phase matrix to 
%   render scatter maps with absolute zero disk I/O and vectorized dilation.
%
%   INPUTS:
%   ps       : pre-loaded PS structure (containing lonlat, etc.)
%   in_ph    : [N_ps x 1] vector of phase or velocity values
%   bg_flag  : 0 (Black background), 1 (White background)
%   col_rg   : Color map limits
%   lon_rg   : Longitude clipping range
%   lat_rg   : Latitude clipping range
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          March 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Zero I/O Memory Architecture: Passes the pre-loaded 'ps' struct directly, 
%      eliminating internal disk loads and drastically reducing execution time.
%   2. Vectorized Rendering: Replaced the O(N) pixel dilation loop with ultra-fast 
%      matrix index expansion and 'accumarray', bypassing MATLAB rendering bottlenecks.
%   3. Mode Pruning: Stripped deprecated DEM, Amplitude, Radar-coord modes, and 
%      external GPS plotting logic to purify the core spatial mapping engine.
%   4. Smart Colormap Routing: Integrated dynamic GMT CPT parsing with intelligent 
%      ocean-clipping (negative value removal) specifically tailored for DEM palettes.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    if nargin < 3, bg_flag = 1; end
    if nargin < 4, col_rg = []; end
    if nargin < 5, lon_rg = []; end
    if nargin < 6, lat_rg = []; end

    % 1. Retrieve global display parameters
    plot_pixel_m = getparm('plot_scatterer_size');
    plot_pixel_size = getparm('plot_pixels_scatterer');
    plot_color_scheme = getparm('plot_color_scheme');
    lonlat_offset = getparm('lonlat_offset');

    % 2. Extract coordinates and apply manual offset
    lonlat = ps.lonlat;
    lonlat(:,1) = lonlat(:,1) + lonlat_offset(1);
    lonlat(:,2) = lonlat(:,2) + lonlat_offset(2);

    % 3. Spatial clipping based on user range
    if ~isempty(lon_rg)
        ix = lonlat(:,1) >= lon_rg(1) & lonlat(:,1) <= lon_rg(2);
        in_ph = in_ph(ix); lonlat = lonlat(ix,:);
    end
    if ~isempty(lat_rg)
        ix = lonlat(:,2) >= lat_rg(1) & lonlat(:,2) <= lat_rg(2);
        in_ph = in_ph(ix); lonlat = lonlat(ix,:);
    end

    if isempty(in_ph)
        ph_lims = [0, 0]; return;
    end

    % 4. Colormap mapping and limits calculation
    if ~isempty(col_rg)
        min_ph = min(col_rg); max_ph = max(col_rg);
    else
        if isreal(in_ph), min_ph = min(in_ph); max_ph = max(in_ph); else, min_ph = 0; max_ph = 0; end
    end
    ph_range = max_ph - min_ph;

    if ~isreal(in_ph) || ph_range == 0
        if ph_range == 0, min_ph = -pi; max_ph = pi; ph_range = 2*pi; end
        in_ph = angle(in_ph); c = hsv(64);
    else
        if strcmpi(plot_color_scheme, 'gray')
            c = flipud(gray(64));  
        elseif strcmpi(plot_color_scheme, 'inflation')
            c = flipud(jet(64)); 
        elseif strcmpi(plot_color_scheme, 'deflation')
            c = jet(64);         
        elseif startsWith(plot_color_scheme, 'GMT_', 'IgnoreCase', true)
            try
                % Smart ocean-clip: remove negative values for DEM palettes
                if ~isempty(regexpi(plot_color_scheme, 'relief|globe|dem|topo|etopo|land|srtm|terra|world|earth|geo'))
                    c = cptcmap(plot_color_scheme, 'ncol', 64, 'clip_neg', true);
                else
                    % Standard load for other GMT palettes (e.g., deformation)
                    c = cptcmap(plot_color_scheme, 'ncol', 64);
                end
            catch
                warning('StaMPS_HPC:CptError', 'Cannot load %s.cpt. Falling back to default jet.', plot_color_scheme);
                c = jet(64);
            end
        else
            c = flipud(jet(64));
        end
    end

    col_ix = double(round(((in_ph - min_ph) * 63 / ph_range) + 1));
    col_ix(col_ix > 64) = 64; col_ix(col_ix < 1) = 1;
    ph_lims = [max_ph, min_ph];

    % Set background color (inject into colormap index 1)
    if bg_flag == 0, c = [[0 0 0]; c]; else, c = [[1 1 1]; c]; end

    % ====================================================================
    % 5. Vectorized Pixel Dilation Engine
    % ====================================================================
    cal0 = mean(lonlat)'; 
    xy_ratio = llh2local([cal0; 0], [cal0+[0.1;0.1]; 0]);
    aspect_ratio = xy_ratio(1) / xy_ratio(2);

    y_posting = plot_pixel_m / plot_pixel_size * 9e-6; 
    x_posting = y_posting / aspect_ratio;

    x = min(lonlat(:,1))-2*x_posting : x_posting : max(lonlat(:,1))+2*x_posting;
    y = min(lonlat(:,2))-2*y_posting : y_posting : max(lonlat(:,2))+2*y_posting;
    Nx = length(x); Ny = length(y);

    % Calculate base grid indices
    demx = round((lonlat(:,1) - x(1)) / x_posting) + 1;
    demy = round((lonlat(:,2) - y(1)) / y_posting) + 1;

    % Remove NaNs and merge overlapping points
    nnix = ~isnan(col_ix);
    demx = demx(nnix); demy = demy(nnix); col_ix = col_ix(nnix);
    [~, uix] = unique([demy, demx], 'rows');
    demx = demx(uix); demy = demy(uix); col_ix = col_ix(uix);

    % Vectorized pixel inflation using accumarray
    if plot_pixel_size == 1
        % Base size: direct mapping
        R = accumarray([demy, demx], col_ix + 1, [Ny, Nx], @max, 1);
    else
        % Dilation size: build neighborhood offset matrix
        margin1 = floor((plot_pixel_size-1)/2);
        margin2 = ceil((plot_pixel_size-1)/2);
        [dC, dR] = meshgrid(-margin1:margin2, -margin1:margin2);
        
        % Matrix broadcast: expand 1 point to N x N pixels
        demy_exp = demy + dR(:)';
        demx_exp = demx + dC(:)';
        col_exp  = repmat(col_ix + 1, 1, numel(dR));
        
        % Flatten arrays
        demy_flat = demy_exp(:); demx_flat = demx_exp(:); col_flat = col_exp(:);
        
        % Boundary safety filter
        valid = demy_flat >= 1 & demy_flat <= Ny & demx_flat >= 1 & demx_flat <= Nx;
        
        % Fast matrix assembly
        R = accumarray([demy_flat(valid), demx_flat(valid)], col_flat(valid), [Ny, Nx], @max, 1);
    end

    % 6. Render and align viewport
    image(x, y, R);
    set(gca, 'YDir', 'normal'); 
    axis tight;
    
    % Set physical aspect ratio
    lon_range = x_posting * Nx; lat_range = y_posting * Ny;
    xy_ratio = llh2local([x(1)+lon_range; y(1)+lat_range; 0], [x(1); y(1); 0]);
    set(gca, 'plotboxaspectratio', [xy_ratio(1)/xy_ratio(2), 1, 1]);
    
    colormap(c);
end