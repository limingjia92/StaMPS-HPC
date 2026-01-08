function plot_patch_kml_gamma(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
% PLOT_PATCH_KML_GAMMA  Generate KML to visualize StaMPS patches from GAMMA source
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 (HPC-Hybrid)
%   License:       GPL v3.0 (Inherited from StaMPS)

% Usage:
%   plot_patch_kml_gamma(data_dir, rg_patches, az_patches)
%   plot_patch_kml_gamma(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
%
% Arguments:
%   data_dir   : Path to the GAMMS source data directory
%   rg_patches : Number of patches in Range
%   az_patches : Number of patches in Azimuth
%   rg_overlap : (Optional) Overlap in Range (Default: 400)
%   az_overlap : (Optional) Overlap in Azimuth (Default: 400)

    % --- 0. 参数默认值处理 ---
    if nargin < 4
        rg_overlap = 400;
    end
    if nargin < 5
        az_overlap = 400;
    end

    fprintf('------------------------------------------------\n');
    fprintf('Generating Patch KML Visualization (GAMMA Mode)\n');
    fprintf('Patches: %d (Rg) x %d (Az)\n', rg_patches, az_patches);
    fprintf('Overlap: %d (Rg), %d (Az)\n', rg_overlap, az_overlap);

    % --- 1. 搜索元数据文件 (Width/Length) ---
    rsc_file = '';
    file_type = ''; % 'par' or 'diff_par'

    % 1.1 尝试搜索 PS 模式路径 (data_dir/rslc/*.par)
    par_files_ps = dir(fullfile(data_dir, 'rslc', '*.par'));
    if ~isempty(par_files_ps)
        rsc_file = fullfile(par_files_ps(1).folder, par_files_ps(1).name);
        file_type = 'par';
        fprintf('Found PS parameter file: %s\n', par_files_ps(1).name);
    end

    % 1.2 如果没找到，尝试搜索 SBAS 模式路径 (data_dir/SMALL_BASELINES/rslc/*.par)
    if isempty(rsc_file)
        par_files_sbas = dir(fullfile(data_dir, 'SMALL_BASELINES', 'rslc', '*.par'));
        if ~isempty(par_files_sbas)
            rsc_file = fullfile(par_files_sbas(1).folder, par_files_sbas(1).name);
            file_type = 'par';
            fprintf('Found SBAS parameter file: %s\n', par_files_sbas(1).name);
        end
    end

    % 1.3 如果还没找到，尝试搜索 GEO 路径 (data_dir/geo/*.diff_par)
    if isempty(rsc_file)
        diff_par_files = dir(fullfile(data_dir, 'geo', '*.diff_par'));
        if ~isempty(diff_par_files)
            rsc_file = fullfile(diff_par_files(1).folder, diff_par_files(1).name);
            file_type = 'diff_par';
            fprintf('Found Diff parameter file: %s\n', diff_par_files(1).name);
        end
    end

    % 错误处理
    if isempty(rsc_file)
        error('Error: Could not find metadata file. Checked:\n - %s\n - %s\n - %s', ...
            fullfile(data_dir, 'rslc', '*.par'), ...
            fullfile(data_dir, 'SMALL_BASELINES', 'rslc', '*.par'), ...
            fullfile(data_dir, 'geo', '*.diff_par'));
    end

    % --- 2. 读取宽长信息 ---
    [img_width, img_length] = read_metadata_dimensions(rsc_file, file_type);
    fprintf('Image Dimensions: %d (Width) x %d (Length)\n', img_width, img_length);

    % --- 3. 自动定位经纬度文件 ---
    geo_dir = fullfile(data_dir, 'geo');
    lon_list = dir(fullfile(geo_dir, '*.lon'));
    lat_list = dir(fullfile(geo_dir, '*.lat'));

    if isempty(lon_list) || isempty(lat_list)
        error('Error: .lon or .lat files not found in %s', geo_dir);
    end

    lon_file = fullfile(lon_list(1).folder, lon_list(1).name);
    lat_file = fullfile(lat_list(1).folder, lat_list(1).name);

    % --- 4. 初始化 KML ---
    kml_filename = 'patches_GAMMA.kml';
    fid_kml = fopen(kml_filename, 'w');
    if fid_kml == -1
        error('Cannot create KML file.');
    end

    write_kml_header(fid_kml);

    % 打开经纬度文件 (使用 Big Endian 读取 float)
    fid_lon = fopen(lon_file, 'r', 'b'); 
    fid_lat = fopen(lat_file, 'r', 'b'); 

    % --- 5. 循环计算 Patch 范围 ---
    width_p = floor(img_width / rg_patches);
    length_p = floor(img_length / az_patches);
    
    patch_count = 0;

    for irg = 1:rg_patches
        for iaz = 1:az_patches
            patch_count = patch_count + 1;

            % --- Calculate Indices (Bash Logic) ---
            start_rg1 = width_p * (irg - 1) + 1;
            end_rg1   = width_p * irg;
            
            start_az1 = length_p * (iaz - 1) + 1;
            end_az1   = length_p * iaz;

            % Apply Overlap
            start_rg = start_rg1 - rg_overlap;
            if start_rg < 1, start_rg = 1; end

            end_rg = end_rg1 + rg_overlap;
            if end_rg > img_width, end_rg = img_width; end

            start_az = start_az1 - az_overlap;
            if start_az < 1, start_az = 1; end

            end_az = end_az1 + az_overlap;
            if end_az > img_length, end_az = img_length; end

            % --- Extract Coordinates for Corners ---
            % Corner order: Top-Left -> Top-Right -> Bottom-Right -> Bottom-Left
            rg_indices = [start_rg, end_rg, end_rg, start_rg];
            az_indices = [start_az, start_az, end_az, end_az];
            
            lons = zeros(1, 4);
            lats = zeros(1, 4);

            for k = 1:4
                lons(k) = get_val_at_pixel(fid_lon, img_width, rg_indices(k), az_indices(k));
                lats(k) = get_val_at_pixel(fid_lat, img_width, rg_indices(k), az_indices(k));
            end
            
            % --- Write to KML ---
            write_kml_placemark(fid_kml, patch_count, lons, lats);
        end
    end

    % --- 6. 清理 ---
    fclose(fid_lon);
    fclose(fid_lat);
    write_kml_footer(fid_kml);
    fclose(fid_kml);

    fprintf('Success! KML file saved to: %s\n', fullfile(pwd, kml_filename));
end

%% --- Helper Functions ---

function val = get_val_at_pixel(fid, width, r, a)
    % 快速读取指定行列的数值 (避免读取整个文件)
    % offset = ((Line-1) * Width + (Col-1)) * 4 bytes
    offset = ((a - 1) * width + (r - 1)) * 4;
    status = fseek(fid, offset, 'bof');
    if status == -1
        warning('Seek failed for pixel (%d, %d)', r, a);
        val = 0;
        return;
    end
    val = fread(fid, 1, 'float', 'b');
end

function [width, len] = read_metadata_dimensions(filepath, file_type)
    % 解析 GAMMA 参数文件获取尺寸
    fid = fopen(filepath, 'r');
    content = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = content{1};
    
    width = 0;
    len = 0;
    
    % 定义要查找的关键字
    if strcmp(file_type, 'par')
        key_width = 'range_samples';
        key_len   = 'azimuth_lines';
    else % diff_par
        key_width = 'map_width';
        key_len   = 'map_azimuth_lines';
    end

    % 解析逻辑
    for i = 1:length(lines)
        line = lines{i};
        if contains(line, key_width)
            % 提取数字: 查找行中的所有数字，取第一个非空的
            tokens = regexp(line, '\d+', 'match');
            if ~isempty(tokens)
                width = str2double(tokens{1});
            end
        end
        if contains(line, key_len)
            tokens = regexp(line, '\d+', 'match');
            if ~isempty(tokens)
                len = str2double(tokens{1});
            end
        end
    end
    
    if width == 0 || len == 0
         error('Could not parse dimensions from %s using keys [%s, %s]', filepath, key_width, key_len);
    end
end

function write_kml_header(fid)
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '<Document>\n');
    fprintf(fid, '  <name>StaMPS Patches</name>\n');
    fprintf(fid, '  <Style id="polygonStyle">\n');
    fprintf(fid, '    <LineStyle>\n');
    fprintf(fid, '      <color>ff0000ff</color>\n'); % Red
    fprintf(fid, '      <width>2</width>\n');
    fprintf(fid, '    </LineStyle>\n');
    fprintf(fid, '    <PolyStyle>\n');
    fprintf(fid, '      <color>00ffffff</color>\n'); % Transparent
    fprintf(fid, '    </PolyStyle>\n');
    fprintf(fid, '  </Style>\n');
end

function write_kml_placemark(fid, idx, lons, lats)
    fprintf(fid, '  <Placemark>\n');
    fprintf(fid, '    <name>PATCH_%d</name>\n', idx);
    fprintf(fid, '    <styleUrl>#polygonStyle</styleUrl>\n');
    fprintf(fid, '    <Polygon>\n');
    fprintf(fid, '      <tessellate>1</tessellate>\n');
    fprintf(fid, '      <outerBoundaryIs>\n');
    fprintf(fid, '        <LinearRing>\n');
    fprintf(fid, '          <coordinates>\n');
    
    % Write coordinates (Lon,Lat,0)
    for i = 1:4
        fprintf(fid, '            %.6f,%.6f,0\n', lons(i), lats(i));
    end
    % Close the loop
    fprintf(fid, '            %.6f,%.6f,0\n', lons(1), lats(1));
    
    fprintf(fid, '          </coordinates>\n');
    fprintf(fid, '        </LinearRing>\n');
    fprintf(fid, '      </outerBoundaryIs>\n');
    fprintf(fid, '    </Polygon>\n');
    fprintf(fid, '  </Placemark>\n');
end

function write_kml_footer(fid)
    fprintf(fid, '</Document>\n');
    fprintf(fid, '</kml>\n');
end