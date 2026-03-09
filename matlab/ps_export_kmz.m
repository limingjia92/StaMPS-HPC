function ps_export_kmz(filename, lonlat, data, varargin)
% PS_EXPORT_KMZ  Export large-scale PS point cloud to a rasterized KMZ file
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
% Usage:
%   ps_export_kmz(filename, lonlat, data)
%   ps_export_kmz(filename, lonlat, data, 'PropertyName', PropertyValue, ...)
%
% Arguments:
%   filename   : Output KMZ filename (e.g., 'insar_result.kmz')
%   lonlat     : N x 2 matrix containing [longitude, latitude]
%   data       : N x 1 vector containing data values for each point
%
% Optional Arguments:
%   'clims'      : Color limits [min_val, max_val] (Default: [min, max])
%   'colormap'   : Colormap string or N x 3 matrix (Default: flipud(jet(64)))
%   'opacity'    : Overlay opacity [0 to 1] (Default: 0.8)
%   'resolution' : Raster spatial resolution in degrees (Default: 0.005)
%   'max_gap'    : Max interpolation gap in pixels (Default: 1)

    %% 1. Parse optional arguments
    clims = [min(data), max(data)];
    cmap = flipud(jet(64)); 
    opacity = 0.8;
    res = 0.005; 
    max_gap = 1;  

    if nargin > 3
        for i = 1:2:length(varargin)
            switch lower(varargin{i})
                case 'clims'
                    clims = varargin{i+1};
                case 'colormap'
                    cmap_input = varargin{i+1};
                    if ischar(cmap_input)
                        cmap = flipud(colormap(cmap_input));
                    else
                        cmap = flipud(cmap_input); 
                    end
                case 'opacity'
                    opacity = varargin{i+1};
                case 'resolution'
                    res = varargin{i+1};
                case 'max_gap'
                    max_gap = varargin{i+1};
                    if ischar(max_gap), max_gap = str2double(max_gap); end
            end
        end
    end

    disp('Generating raster grid and adaptive mask...');
    
    %% 2. Create spatial grid and adaptive mask
    min_lon = min(lonlat(:,1)); max_lon = max(lonlat(:,1));
    min_lat = min(lonlat(:,2)); max_lat = max(lonlat(:,2));

    lon_grid = min_lon : res : max_lon;
    lat_grid = min_lat : res : max_lat;
    [X, Y] = meshgrid(lon_grid, lat_grid);

    % Build validity mask
    col_idx = round((lonlat(:,1) - min_lon) / res) + 1;
    row_idx = round((lonlat(:,2) - min_lat) / res) + 1;
    col_idx = max(1, min(length(lon_grid), col_idx));
    row_idx = max(1, min(length(lat_grid), row_idx));
    
    point_mask = false(length(lat_grid), length(lon_grid));
    linear_idx = sub2ind(size(point_mask), row_idx, col_idx);
    point_mask(linear_idx) = true;
    
    % Morphological dilation to fill small gaps
    [xx, yy] = meshgrid(-max_gap:max_gap, -max_gap:max_gap);
    kernel = double(xx.^2 + yy.^2 <= max_gap^2);
    smoothed_mask = conv2(double(point_mask), kernel, 'same');
    valid_mask = smoothed_mask > 0;
    
    %% 3. Spatial interpolation
    disp('Performing spatial interpolation...');
    F = scatteredInterpolant(lonlat(:,1), lonlat(:,2), double(data), 'linear', 'nearest');
    Z = F(X, Y);
    
    % Crop invalid areas (e.g., oceans)
    Z(~valid_mask) = NaN; 

    %% 4. Map data to RGB colormap
    disp('Rendering image and alpha channel...');
    Z_norm = (Z - clims(1)) / (clims(2) - clims(1));
    Z_norm(Z_norm < 0) = 0;
    Z_norm(Z_norm > 1) = 1;

    nColors = size(cmap, 1);
    Z_idx = round(Z_norm * (nColors - 1)) + 1;

    RGB = zeros(size(Z, 1), size(Z, 2), 3);
    for i = 1:3
        color_channel = cmap(:, i);
        temp_img = zeros(size(Z));
        valid_pixels = ~isnan(Z);
        temp_img(valid_pixels) = color_channel(Z_idx(valid_pixels));
        RGB(:, :, i) = temp_img;
    end

    % Flip up-down to match image coordinate system
    RGB = flipud(RGB);
    Z = flipud(Z);

    % Construct Alpha channel
    Alpha = opacity * ones(size(Z));
    Alpha(isnan(Z)) = 0;

    %% 5. Convert to 8-bit to prevent Google Earth render bugs
    RGB_uint8 = uint8(RGB * 255);
    Alpha_uint8 = uint8(Alpha * 255);

    %% 6. Generate temporary files and pack into KMZ
    disp('Packing KMZ file...');
    png_filename = 'overlay_temp.png';
    kml_filename = 'doc.kml';
    
    imwrite(RGB_uint8, png_filename, 'Alpha', Alpha_uint8);

    fid = fopen(kml_filename, 'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '  <Folder>\n');
    fprintf(fid, '    <name>InSAR Raster Result</name>\n');
    fprintf(fid, '    <GroundOverlay>\n');
    fprintf(fid, '      <name>Deformation Rate</name>\n');
    fprintf(fid, '      <color>ffffffff</color>\n');
    fprintf(fid, '      <Icon>\n');
    fprintf(fid, '        <href>%s</href>\n', png_filename);
    fprintf(fid, '      </Icon>\n');
    fprintf(fid, '      <LatLonBox>\n');
    fprintf(fid, '        <north>%.6f</north>\n', max_lat);
    fprintf(fid, '        <south>%.6f</south>\n', min_lat);
    fprintf(fid, '        <east>%.6f</east>\n', max_lon);
    fprintf(fid, '        <west>%.6f</west>\n', min_lon);
    fprintf(fid, '      </LatLonBox>\n');
    fprintf(fid, '    </GroundOverlay>\n');
    fprintf(fid, '  </Folder>\n');
    fprintf(fid, '</kml>\n');
    fclose(fid);

    zip('temp_archive.zip', {png_filename, kml_filename});
    movefile('temp_archive.zip', filename, 'f');
    delete(png_filename);
    delete(kml_filename);

    disp(['Done! File saved to: ', filename]);
end