function plot_patch_kml_isce(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
% PLOT_PATCH_KML_ISCE  Generate KML and GMT text to visualize Patches (ISCE)
%
%   Usage:
%       plot_patch_kml_isce(data_dir, rg_patches, az_patches)
%       plot_patch_kml_isce(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
%
%   Inputs:
%       data_dir   - Directory containing the ISCE processed data.
%                    Expected Data Directory Structure must include:
%                    * Dimension files: 'width.txt' AND 'len.txt'
%                    * Geometry files: 'lon.raw' AND 'lat.raw'
%       rg_patches - Number of patches in the range direction.
%       az_patches - Number of patches in the azimuth direction.
%       rg_overlap - (Optional) Number of overlapping pixels between adjacent range patches (Default: 400).
%       az_overlap - (Optional) Number of overlapping pixels between adjacent azimuth patches (Default: 400).
%
%   Outputs:
%       Generates the following files in the current working directory:
%       - patches_ISCE.kml : KML file for visualization in Google Earth or other GIS software.
%       - patches_ISCE.txt : Multi-segment polygon text file for GMT plotting.
%
%   Example:
%       % Generate visualization for 2 range patches and 5 azimuth patches with default 400px overlap
%       plot_patch_kml_isce('/path/to/data_dir', 2, 5);
%
%       % Specify custom overlap of 300 pixels in both directions
%       plot_patch_kml_isce('/path/to/data_dir', 2, 5, 300, 300);
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       2.0 (Added GMT output & KML LookAt)
%   ======================================================================

    if nargin < 4, rg_overlap = 400; end
    if nargin < 5, az_overlap = 400; end

    fprintf('------------------------------------------------\n');
    fprintf('Generating Patch Visualization (ISCE Mode)\n');
    fprintf('Data Dir: %s\n', data_dir);
    fprintf('Patches : %d (Rg) x %d (Az)\n', rg_patches, az_patches);

    % --- 1. Read Width and Length Information ---
    width_file = fullfile(data_dir, 'width.txt');
    len_file   = fullfile(data_dir, 'len.txt');

    if ~isfile(width_file) || ~isfile(len_file)
        error('Error: width.txt or len.txt not found in %s.', data_dir);
    end
    img_width  = load(width_file);
    img_length = load(len_file);

    % --- 2. Locate Longitude/Latitude Files ---
    lon_file = fullfile(data_dir, 'lon.raw');
    lat_file = fullfile(data_dir, 'lat.raw');

    if ~isfile(lon_file) || ~isfile(lat_file)
        error('Error: lon.raw or lat.raw not found.');
    end

    fid_lon = fopen(lon_file, 'r', 'l'); 
    fid_lat = fopen(lat_file, 'r', 'l'); 

    % --- 3. Pre-calculate Global Bounding Box for KML <LookAt> ---
    global_corners = estimate_corners_affine(fid_lon, fid_lat, img_width, ...
                                             1, img_width, 1, img_length, ...
                                             1, img_width, 1, img_length);
    if ~isempty(global_corners)
        min_lon = min(global_corners(:,1)); max_lon = max(global_corners(:,1));
        min_lat = min(global_corners(:,2)); max_lat = max(global_corners(:,2));
        center_lon = (min_lon + max_lon) / 2;
        center_lat = (min_lat + max_lat) / 2;
        % calculate look at angle
        lookat_range = max(max_lon - min_lon, max_lat - min_lat) * 111000 * 1.5; 
    else
        center_lon = 0; center_lat = 0; lookat_range = 100000;
    end

    % --- 4. Initialize Output Files (KML + GMT Text) ---
    kml_filename = 'patches_ISCE.kml';
    txt_filename = 'patches_ISCE.txt';
    
    fid_kml = fopen(kml_filename, 'w');
    fid_txt = fopen(txt_filename, 'w');
    
    write_kml_header(fid_kml, center_lon, center_lat, lookat_range);

    % --- 5. Loop to Calculate Patches ---
    width_p = floor(img_width / rg_patches);
    length_p = floor(img_length / az_patches);
    patch_count = 0;

    for irg = 1:rg_patches
        for iaz = 1:az_patches
            patch_count = patch_count + 1;

            % Grid based indices
            start_rg_noover = width_p * (irg - 1) + 1;
            end_rg_noover   = width_p * irg;
            start_az_noover = length_p * (iaz - 1) + 1;
            end_az_noover   = length_p * iaz;

            % Apply Overlap (Theoretical limits)
            start_rg = start_rg_noover - rg_overlap;
            end_rg   = end_rg_noover + rg_overlap;
            start_az = start_az_noover - az_overlap;
            end_az   = end_az_noover + az_overlap;

            % Clamp to image limits for data reading
            read_start_rg = max(1, start_rg); 
            read_end_rg   = min(img_width, end_rg);
            read_start_az = max(1, start_az); 
            read_end_az   = min(img_length, end_az);

            % Affine fitting to estimate corners
            corners = estimate_corners_affine(fid_lon, fid_lat, img_width, ...
                                              read_start_rg, read_end_rg, ...
                                              read_start_az, read_end_az, ...
                                              start_rg, end_rg, start_az, end_az);
            
            % Write to outputs
            if ~isempty(corners)
                lons = corners(:, 1);
                lats = corners(:, 2);
                
                % write KML
                write_kml_placemark(fid_kml, patch_count, lons, lats);
                
                % write GMT data file
                fprintf(fid_txt, '> PATCH_%d\n', patch_count);
                for i = 1:4
                    fprintf(fid_txt, '%.6f %.6f\n', lons(i), lats(i));
                end
                % polygon
                fprintf(fid_txt, '%.6f %.6f\n', lons(1), lats(1));
                
                fprintf('Patch %d: Processed.\n', patch_count);
            else
                fprintf('Warning: Patch %d has insufficient valid data (skipped).\n', patch_count);
            end
        end
    end

    % --- 6. Cleanup ---
    fclose(fid_lon); fclose(fid_lat);
    write_kml_footer(fid_kml); fclose(fid_kml);
    fclose(fid_txt);
    fprintf('Success! Saved to:\n - %s\n - %s\n', kml_filename, txt_filename);
end

%% --- Helper Functions ---
function predicted_corners = estimate_corners_affine(fid_lon, fid_lat, full_width, ...
                                                     r_min, r_max, a_min, a_max, ...
                                                     target_r_min, target_r_max, target_a_min, target_a_max)
    predicted_corners = [];
    n_samples = 20; 
    step_r = max(1, floor((r_max - r_min) / n_samples));
    step_a = max(1, floor((a_max - a_min) / n_samples));
    r_coords = []; a_coords = []; val_lons = []; val_lats = [];
    
    for a = a_min:step_a:a_max
        for r = r_min:step_r:r_max
             l_val = get_val_at_pixel(fid_lon, full_width, r, a);
             L_val = get_val_at_pixel(fid_lat, full_width, r, a);
             if abs(l_val) > 1e-6 && abs(L_val) > 1e-6
                 r_coords = [r_coords; r]; a_coords = [a_coords; a]; %#ok<AGROW>
                 val_lons = [val_lons; l_val]; val_lats = [val_lats; L_val]; %#ok<AGROW>
             end
        end
    end
    
    if length(r_coords) < 10, return; end
    
    A = [r_coords, a_coords, ones(size(r_coords))];
    coeffs_lon = A \ val_lons; 
    coeffs_lat = A \ val_lats;
    
    target_corners_r = [target_r_min; target_r_max; target_r_max; target_r_min];
    target_corners_a = [target_a_min; target_a_min; target_a_max; target_a_max];
    DesignMatrix = [target_corners_r, target_corners_a, ones(4,1)];
    
    pred_lons = DesignMatrix * coeffs_lon;
    pred_lats = DesignMatrix * coeffs_lat;
    predicted_corners = [pred_lons, pred_lats];
end

function val = get_val_at_pixel(fid, width, r, a)
    offset = ((a - 1) * width + (r - 1)) * 4;
    fseek(fid, offset, 'bof');
    val = fread(fid, 1, 'float', 'l');
    if isempty(val), val = 0; end
end

function write_kml_header(fid, c_lon, c_lat, range)
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '<Document>\n');
    fprintf(fid, '  <name>StaMPS Patches (ISCE)</name>\n');
    
    fprintf(fid, '  <LookAt>\n');
    fprintf(fid, '    <longitude>%.6f</longitude>\n', c_lon);
    fprintf(fid, '    <latitude>%.6f</latitude>\n', c_lat);
    fprintf(fid, '    <altitude>0</altitude>\n');
    fprintf(fid, '    <heading>0</heading>\n');
    fprintf(fid, '    <tilt>0</tilt>\n');
    fprintf(fid, '    <range>%.0f</range>\n', range);
    fprintf(fid, '  </LookAt>\n');
    
    fprintf(fid, '  <Style id="centeredLabelStyle">\n');
    fprintf(fid, '    <IconStyle><scale>0</scale></IconStyle>\n');
    fprintf(fid, '    <LabelStyle>\n');
    fprintf(fid, '      <scale>1.5</scale>\n');
    fprintf(fid, '      <color>ff0000ff</color>\n'); 
    fprintf(fid, '    </LabelStyle>\n');
    fprintf(fid, '    <LineStyle>\n');
    fprintf(fid, '      <color>ff0000ff</color>\n'); 
    fprintf(fid, '      <width>2</width>\n');  
    fprintf(fid, '    </LineStyle>\n');
    fprintf(fid, '    <PolyStyle><fill>0</fill></PolyStyle>\n');
    fprintf(fid, '  </Style>\n');
end

function write_kml_placemark(fid, idx, lons, lats)
    lon_center = mean(lons); lat_center = mean(lats);
    fprintf(fid, '  <Placemark>\n');
    fprintf(fid, '    <name>PATCH_%d</name>\n', idx);
    fprintf(fid, '    <styleUrl>#centeredLabelStyle</styleUrl>\n'); 
    fprintf(fid, '    <MultiGeometry>\n');
    fprintf(fid, '      <Polygon>\n');
    fprintf(fid, '        <tessellate>1</tessellate>\n');
    fprintf(fid, '        <outerBoundaryIs>\n');
    fprintf(fid, '          <LinearRing>\n');
    fprintf(fid, '            <coordinates>\n');
    for i = 1:4
        fprintf(fid, '              %.6f,%.6f,0\n', lons(i), lats(i));
    end
    fprintf(fid, '              %.6f,%.6f,0\n', lons(1), lats(1));
    fprintf(fid, '            </coordinates>\n');
    fprintf(fid, '          </LinearRing>\n');
    fprintf(fid, '        </outerBoundaryIs>\n');
    fprintf(fid, '      </Polygon>\n');
    fprintf(fid, '      <Point>\n');
    fprintf(fid, '        <coordinates>%.6f,%.6f,0</coordinates>\n', lon_center, lat_center);
    fprintf(fid, '      </Point>\n');
    fprintf(fid, '    </MultiGeometry>\n');
    fprintf(fid, '  </Placemark>\n');
end

function write_kml_footer(fid)
    fprintf(fid, '</Document>\n</kml>\n');
end