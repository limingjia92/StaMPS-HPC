function plot_patch_kml_gamma(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
% PLOT_PATCH_KML_GAMMA  Generate KML and GMT text to visualize Patches (GAMMA)
%
%   Usage:
%       plot_patch_kml_gamma(data_dir, rg_patches, az_patches)
%       plot_patch_kml_gamma(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
%
%   Inputs:
%       data_dir   - Directory containing the GAMMA processed data. 
%                    Expected Data Directory Structure must include:
%                    * Metadata file: 'rslc/*.par' OR 'SMALL_BASELINES/rslc/*.par' OR 'geo/*.diff_par'
%                    * Geometry files: 'geo/*.lon' AND 'geo/*.lat'
%       rg_patches - Number of patches in the range direction.
%       az_patches - Number of patches in the azimuth direction.
%       rg_overlap - (Optional) Number of overlapping pixels between adjacent range patches (Default: 400).
%       az_overlap - (Optional) Number of overlapping pixels between adjacent azimuth patches (Default: 400).
%
%   Outputs:
%       Generates the following files in the current working directory:
%       - patches_GAMMA.kml : KML file for visualization in Google Earth or other GIS software.
%       - patches_GAMMA.txt : Multi-segment polygon text file for GMT plotting.
%
%   Example:
%       % Generate visualization for 2 range patches and 5 azimuth patches with default 400px overlap
%       plot_patch_kml_gamma('/path/to/data_dir', 2, 5);
%
%       % Specify custom overlap of 300 pixels in both directions
%       plot_patch_kml_gamma('/path/to/data_dir', 2, 5, 300, 300);
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
    fprintf('Generating Patch Visualization (GAMMA Mode)\n');
    fprintf('Patches: %d (Rg) x %d (Az)\n', rg_patches, az_patches);
    
    % --- 1. Search Metadata File ---
    rsc_file = ''; file_type = '';
    search_paths = {fullfile(data_dir, 'rslc', '*.par'), ...
                    fullfile(data_dir, 'SMALL_BASELINES', 'rslc', '*.par'), ...
                    fullfile(data_dir, 'geo', '*.diff_par')};
    
    for i = 1:length(search_paths)
        files = dir(search_paths{i});
        if ~isempty(files)
            rsc_file = fullfile(files(1).folder, files(1).name);
            file_type = 'par'; if i == 3, file_type = 'diff_par'; end
            break;
        end
    end
    if isempty(rsc_file), error('Error: Could not find metadata file.'); end

    % --- 2. Read Dimensions ---
    [img_width, img_length] = read_metadata_dimensions(rsc_file, file_type);

    % --- 3. Locate Longitude/Latitude Files ---
    geo_dir = fullfile(data_dir, 'geo');
    lon_list = dir(fullfile(geo_dir, '*.lon'));
    lat_list = dir(fullfile(geo_dir, '*.lat'));
    if isempty(lon_list) || isempty(lat_list), error('.lon or .lat files not found.'); end

    lon_file = fullfile(lon_list(1).folder, lon_list(1).name);
    lat_file = fullfile(lat_list(1).folder, lat_list(1).name);

    fid_lon = fopen(lon_file, 'r', 'b'); 
    fid_lat = fopen(lat_file, 'r', 'b'); 

    % --- 4. Pre-calculate Global Bounding Box for KML <LookAt> ---
    % Test the 4 extreme corners of the image
    c_lons = zeros(1,4); c_lats = zeros(1,4);
    r_ext = [1, img_width, img_width, 1];
    a_ext = [1, 1, img_length, img_length];
    for k = 1:4
        c_lons(k) = get_val_at_pixel(fid_lon, img_width, r_ext(k), a_ext(k));
        c_lats(k) = get_val_at_pixel(fid_lat, img_width, r_ext(k), a_ext(k));
    end
    
    valid = (abs(c_lons) > 1e-6);
    if any(valid)
        min_lon = min(c_lons(valid)); max_lon = max(c_lons(valid));
        min_lat = min(c_lats(valid)); max_lat = max(c_lats(valid));
        center_lon = (min_lon + max_lon) / 2;
        center_lat = (min_lat + max_lat) / 2;
        lookat_range = max(max_lon - min_lon, max_lat - min_lat) * 111000 * 1.5;
    else
        center_lon = 0; center_lat = 0; lookat_range = 100000;
    end

    % --- 5. Initialize Output Files ---
    kml_filename = 'patches_GAMMA.kml';
    txt_filename = 'patches_GAMMA.txt';
    fid_kml = fopen(kml_filename, 'w');
    fid_txt = fopen(txt_filename, 'w');

    write_kml_header(fid_kml, center_lon, center_lat, lookat_range);

    % --- 6. Loop to Calculate Patches ---
    width_p = floor(img_width / rg_patches);
    length_p = floor(img_length / az_patches);
    patch_count = 0;

    for irg = 1:rg_patches
        for iaz = 1:az_patches
            patch_count = patch_count + 1;

            start_rg1 = width_p * (irg - 1) + 1;
            end_rg1   = width_p * irg;
            start_az1 = length_p * (iaz - 1) + 1;
            end_az1   = length_p * iaz;

            % Apply Overlap
            start_rg = max(1, start_rg1 - rg_overlap);
            end_rg   = min(img_width, end_rg1 + rg_overlap);
            start_az = max(1, start_az1 - az_overlap);
            end_az   = min(img_length, end_az1 + az_overlap);

            % Corner order: TL -> TR -> BR -> BL
            rg_indices = [start_rg, end_rg, end_rg, start_rg];
            az_indices = [start_az, start_az, end_az, end_az];
            
            lons = zeros(1, 4); lats = zeros(1, 4);
            for k = 1:4
                lons(k) = get_val_at_pixel(fid_lon, img_width, rg_indices(k), az_indices(k));
                lats(k) = get_val_at_pixel(fid_lat, img_width, rg_indices(k), az_indices(k));
            end
            
            % Write to KML
            write_kml_placemark(fid_kml, patch_count, lons, lats);
            
            % Write to GMT Multi-segment text file
            fprintf(fid_txt, '> PATCH_%d\n', patch_count);
            for i = 1:4
                fprintf(fid_txt, '%.6f %.6f\n', lons(i), lats(i));
            end
            fprintf(fid_txt, '%.6f %.6f\n', lons(1), lats(1)); % Close poly
        end
    end

    % --- 7. Cleanup ---
    fclose(fid_lon); fclose(fid_lat);
    write_kml_footer(fid_kml); fclose(fid_kml); fclose(fid_txt);

    fprintf('Success! Saved to:\n - %s\n - %s\n', kml_filename, txt_filename);
end

%% --- Helper Functions ---
function val = get_val_at_pixel(fid, width, r, a)
    offset = ((a - 1) * width + (r - 1)) * 4;
    status = fseek(fid, offset, 'bof');
    if status == -1, val = 0; return; end
    val = fread(fid, 1, 'float', 'b');
end

function [width, len] = read_metadata_dimensions(filepath, file_type)
    fid = fopen(filepath, 'r');
    content = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
    lines = content{1}; width = 0; len = 0;
    
    if strcmp(file_type, 'par')
        key_width = 'range_samples'; key_len = 'azimuth_lines';
    else
        key_width = 'map_width'; key_len = 'map_azimuth_lines';
    end

    for i = 1:length(lines)
        if contains(lines{i}, key_width)
            t = regexp(lines{i}, '\d+', 'match'); if ~isempty(t), width = str2double(t{1}); end
        end
        if contains(lines{i}, key_len)
            t = regexp(lines{i}, '\d+', 'match'); if ~isempty(t), len = str2double(t{1}); end
        end
    end
end

function write_kml_header(fid, c_lon, c_lat, range)
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '<Document>\n');
    fprintf(fid, '  <name>StaMPS Patches (GAMMA)</name>\n');
    
    % LookAt 标签
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