function plot_patch_kml_isce(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
% PLOT_PATCH_KML_ISCE  Generate KML to visualize StaMPS patches from ISCE source
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 (HPC-Hybrid)
%   License:       GPL v3.0 (Inherited from StaMPS)

% Usage:
%   plot_patch_kml_isce(data_dir, rg_patches, az_patches)
%   plot_patch_kml_isce(data_dir, rg_patches, az_patches, rg_overlap, az_overlap)
%
% Arguments:
%   data_dir   : Path to the GAMMS source data directory
%   rg_patches : Number of patches in Range
%   az_patches : Number of patches in Azimuth
%   rg_overlap : (Optional) Overlap in Range (Default: 400)
%   az_overlap : (Optional) Overlap in Azimuth (Default: 400)

    % --- 0. 参数默认值处理 ---
    if nargin < 4, rg_overlap = 400; end
    if nargin < 5, az_overlap = 400; end

    fprintf('------------------------------------------------\n');
    fprintf('Generating Patch KML Visualization (ISCE Mode)\n');
    fprintf('Data Dir: %s\n', data_dir);
    fprintf('Patches : %d (Rg) x %d (Az)\n', rg_patches, az_patches);

    % --- 1. 读取宽长信息 ---
    width_file = fullfile(data_dir, 'width.txt');
    len_file   = fullfile(data_dir, 'len.txt');

    if ~isfile(width_file) || ~isfile(len_file)
        error('Error: width.txt or len.txt not found in %s.', data_dir);
    end

    img_width  = load(width_file);
    img_length = load(len_file);

    % --- 2. 定位经纬度文件 ---
    lon_file = fullfile(data_dir, 'lon.raw');
    lat_file = fullfile(data_dir, 'lat.raw');

    if ~isfile(lon_file) || ~isfile(lat_file)
        error('Error: lon.raw or lat.raw not found.');
    end

    % --- 3. 初始化 KML ---
    kml_filename = 'patches_ISCE.kml';
    fid_kml = fopen(kml_filename, 'w');
    write_kml_header(fid_kml);

    % 打开文件 (ISCE Little Endian)
    fid_lon = fopen(lon_file, 'r', 'l'); 
    fid_lat = fopen(lat_file, 'r', 'l'); 

    % --- 4. 循环计算 Patch ---
    width_p = floor(img_width / rg_patches);
    length_p = floor(img_length / az_patches);
    
    patch_count = 0;

    for irg = 1:rg_patches
        for iaz = 1:az_patches
            patch_count = patch_count + 1;

            % --- Calculate Indices (with Overlap) ---
            % Grid based
            start_rg_noover = width_p * (irg - 1) + 1;
            end_rg_noover   = width_p * irg;
            start_az_noover = length_p * (iaz - 1) + 1;
            end_az_noover   = length_p * iaz;

            % Apply Overlap (Theoretical limits)
            start_rg = start_rg_noover - rg_overlap;
            end_rg   = end_rg_noover + rg_overlap;
            start_az = start_az_noover - az_overlap;
            end_az   = end_az_noover + az_overlap;

            % Clamp to image limits for data reading, but keep theoretical shape for geometry?
            % Usually we want to visualize the processing area. Let's clamp for validity.
            read_start_rg = max(1, start_rg); 
            read_end_rg   = min(img_width, end_rg);
            read_start_az = max(1, start_az); 
            read_end_az   = min(img_length, end_az);

            % --- 核心：仿射拟合推算角点 ---
            corners = estimate_corners_affine(fid_lon, fid_lat, img_width, ...
                                              read_start_rg, read_end_rg, ...
                                              read_start_az, read_end_az, ...
                                              start_rg, end_rg, start_az, end_az);
            
            % --- Write to KML ---
            if ~isempty(corners)
                % corners is 4x2 [lon, lat]
                lons = corners(:, 1);
                lats = corners(:, 2);
                
                write_kml_placemark(fid_kml, patch_count, lons, lats);
                fprintf('Patch %d: Processed.\n', patch_count);
            else
                fprintf('Warning: Patch %d has insufficient valid data for fit (skipped).\n', patch_count);
            end
        end
    end

    % --- 5. 清理 ---
    fclose(fid_lon); fclose(fid_lat);
    write_kml_footer(fid_kml); fclose(fid_kml);
    fprintf('Success! Saved to: %s\n', fullfile(pwd, kml_filename));
end

%% --- Helper Functions ---

function predicted_corners = estimate_corners_affine(fid_lon, fid_lat, full_width, ...
                                                     r_min, r_max, a_min, a_max, ...
                                                     target_r_min, target_r_max, target_a_min, target_a_max)
    % ESTIMATE_CORNERS_AFFINE
    % 1. 在 Patch 范围内采样有效点
    % 2. 拟合平面模型 Lon = f(r, a), Lat = g(r, a)
    % 3. 推算目标角点的经纬度
    
    predicted_corners = [];
    
    % Step 1: 稀疏采样 (Grid Sampling)
    % 为了速度，不在每个像素都读。在 Range 和 Azimuth 方向各采 ~20 个点
    n_samples = 20; 
    step_r = max(1, floor((r_max - r_min) / n_samples));
    step_a = max(1, floor((a_max - a_min) / n_samples));
    
    r_coords = [];
    a_coords = [];
    val_lons = [];
    val_lats = [];
    
    % 遍历采样点
    for a = a_min:step_a:a_max
        for r = r_min:step_r:r_max
             l_val = get_val_at_pixel(fid_lon, full_width, r, a);
             L_val = get_val_at_pixel(fid_lat, full_width, r, a);
             
             % 过滤无效值 (0)
             if abs(l_val) > 1e-6 && abs(L_val) > 1e-6
                 r_coords = [r_coords; r]; %#ok<AGROW>
                 a_coords = [a_coords; a]; %#ok<AGROW>
                 val_lons = [val_lons; l_val]; %#ok<AGROW>
                 val_lats = [val_lats; L_val]; %#ok<AGROW>
             end
        end
    end
    
    % 如果有效点太少，无法拟合
    if length(r_coords) < 10
        return;
    end
    
    % Step 2: 最小二乘法拟合 (Least Squares Fit)
    % Model: Val = c1 * r + c2 * a + c3
    % Matrix: [r, a, 1] * [c1; c2; c3] = Val
    
    A = [r_coords, a_coords, ones(size(r_coords))];
    
    % Fit Longitude
    coeffs_lon = A \ val_lons; % Solve linear system
    
    % Fit Latitude
    coeffs_lat = A \ val_lats;
    
    % Step 3: 推算四个角点 (Predict Corners)
    % Order: TL, TR, BR, BL (ensure counter-clockwise or loop)
    % Corner 1: (min_r, min_a) -> Top-Left (usually)
    % Corner 2: (max_r, min_a) -> Top-Right
    % Corner 3: (max_r, max_a) -> Bottom-Right
    % Corner 4: (min_r, max_a) -> Bottom-Left
    
    % 注意：这里使用的是 target_r/a (包含 overlap 的理论边界)
    % 即使这些点在原始数据中是 0 或超出范围，拟合公式也能给出正确的推算值
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

function write_kml_header(fid)
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid, '<Document>\n<name>StaMPS Patches (Skewed)</name>\n');
    fprintf(fid, '<Style id="pStyle"><LineStyle><color>ff0000ff</color><width>2</width></LineStyle><PolyStyle><color>00ffffff</color></PolyStyle></Style>\n');
end

function write_kml_placemark(fid, idx, lons, lats)
    fprintf(fid, '<Placemark><name>P_%d</name><styleUrl>#pStyle</styleUrl><Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>\n', idx);
    % Loop corners
    for i = 1:4
        fprintf(fid, '%.6f,%.6f,0\n', lons(i), lats(i));
    end
    % Close loop
    fprintf(fid, '%.6f,%.6f,0\n', lons(1), lats(1)); 
    fprintf(fid, '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n');
end

function write_kml_footer(fid)
    fprintf(fid, '</Document></kml>\n');
end