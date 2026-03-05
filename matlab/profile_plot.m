function profile_plot(h_fig)
%PROFILE_PLOT Interactive swath profile plotting 
%   Triggered by a UI button dynamically generated in ps_plot.
%
%   INPUTS:
%   h_fig    : Handle of the main scatter plot figure containing bound data.
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
%   1. Interactive UX: Implemented mouse-driven dual-point selection with 
%      dynamic line rendering and global handle tracking via appdata.
%   2. Vectorized Binning: Eradicated the legacy O(N) loop, utilizing 
%      'discretize' and 'accumarray' for instantaneous swath aggregation.
%   3. Dual-Axis Integration: Automatically fetches topography data to render 
%      velocity and elevation synchronously on a modernized dual-Y axis.
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, December 2010
%   ======================================================================


    % 1. Retrieve bound TS data file
    ts_matname = getappdata(h_fig, 'ts_matname');
    if isempty(ts_matname) || ~exist(ts_matname, 'file')
        errordlg('Cannot find bound TS data file. Please regenerate.', 'Data Error');
        return;
    end
    
    % 2. Get user-defined extraction radius (Swath half-width)
    mEditBox = findobj(h_fig, 'Tag', 'RadiusEditBox');
    if isempty(mEditBox)
        radiusfactor = 500;
    else
        radius_vals = str2num(char(get(mEditBox,'String'))); %#ok<ST2NM>
        if isempty(radius_vals), radiusfactor = 500; else, radiusfactor = radius_vals(1); end
    end
    
    % Get user-defined number of bins
    bEditBox = findobj(h_fig, 'Tag', 'BinsEditBox');
    if isempty(bEditBox)
        bins = 20;
    else
        bin_vals = str2num(char(get(bEditBox,'String'))); %#ok<ST2NM>
        % 防呆设计：bins必须是大于等于2的整数，否则无视输入使用默认值20
        if isempty(bin_vals)
            bins = 20; 
        else
            bins = max(2, round(bin_vals(1))); 
        end
    end

    % 3. Load exclusive data for this figure
    ts_data = load(ts_matname);
    lonlat = ts_data.lonlat;
    
    % Check for velocity field result
    if isfield(ts_data, 'ph_all')
        v_all = ts_data.ph_all;
    else
        errordlg('Variable ''ph_all'' (velocity field) not found in bound data.', 'Data Error');
        return;
    end
    
    momfig_name = ['(', get(h_fig, 'Name'), ')'];
    
    % Dynamically check and load Topography (Elevation)
    try
        load('psver.mat', 'psver');
        hgtname = ['./hgt', num2str(psver), '.mat'];
        has_hgt = exist(hgtname, 'file');
        if has_hgt
            hgt_data = load(hgtname);
            h_all = hgt_data.hgt;
        end
    catch
        has_hgt = false;
    end
    
    % 4. Activate main figure and await dual user selections
    figure(h_fig);
    
    % Clear previous residual markers and lines (shared handle space with TS)
    old_handles = getappdata(h_fig, 'ts_marker_handle');
    if ~isempty(old_handles) && all(isgraphics(old_handles))
        delete(old_handles); 
    end

    hold(gca, 'on');
    disp(['Select 1ST point (Start) for Profile (Radius: ', num2str(radiusfactor), ' m)...']);
    try
        [lon1, lat1] = ginput(1);
        h_m1 = plot(gca, lon1, lat1, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        setappdata(h_fig, 'ts_marker_handle', h_m1); % Prevent ghosting
        disp(['  -> Start Point: ', num2str(lon1), ', ', num2str(lat1)]);
    catch
        disp('Point 1 selection cancelled.'); hold(gca, 'off'); return;
    end
    
    disp('Select 2ND point (End) for Profile...');
    try
        [lon2, lat2] = ginput(1);
        h_m2 = plot(gca, lon2, lat2, 'bx', 'MarkerSize', 12, 'LineWidth', 2.5); 
        
        % Draw black connecting line
        h_line = plot(gca, [lon1, lon2], [lat1, lat2], 'k-', 'LineWidth', 2);
        
        % Update appdata with all three handles
        setappdata(h_fig, 'ts_marker_handle', [h_m1, h_m2, h_line]); 
        disp(['  -> End Point: ', num2str(lon2), ', ', num2str(lat2)]);
    catch
        disp('Point 2 selection cancelled.'); hold(gca, 'off'); return;
    end
    hold(gca, 'off');

    % 5. Coordinate Transformation (Local projection in km)
    lonlatXY = llh2local(lonlat', [lon1; lat1]);
    P1 = llh2local([lon1; lat1], [lon1; lat1]); % Origin [0; 0]
    P2 = llh2local([lon2; lat2], [lon1; lat1]);
    
    % Calculate rotation angle to align profile line horizontally
    dx = P2(1) - P1(1); 
    dy = P2(2) - P1(2);
    alpha = atan2(dy, dx);
    
    % Rotation matrix
    rot = [cos(alpha), -sin(alpha); 
           sin(alpha),  cos(alpha)];
           
    XY_new = rot \ lonlatXY;
    P2_new = rot \ P2;
    profile_length = P2_new(1); % Total profile length in km

    % 6. Swath Filtering (Extract PS points within tube)
    R_km = radiusfactor / 1000;
    
    % Filter by perpendicular distance (Y) and bounded by X (with slight buffer)
    ix_tube = abs(XY_new(2,:)) <= R_km & ...
              XY_new(1,:) >= -R_km & ...
              XY_new(1,:) <= profile_length + R_km;
              
    n_pts_used = sum(ix_tube);
    if n_pts_used == 0
        errordlg('No PS points found within the profile swath. Try a larger radius.', 'No Data');
        return;
    end
    disp([num2str(n_pts_used), ' PS points extracted along the profile.']);
    
    X_tube = XY_new(1, ix_tube)';
    v_tube = v_all(ix_tube);
    if has_hgt, h_tube = h_all(ix_tube); end

    % 7. Vectorized Binning (Replaced legacy for-loop with accumarray)
    [bin_idx, edges] = discretize(X_tube, bins);
    bin_centers = edges(1:end-1) + diff(edges)/2;
    
    valid_idx = ~isnan(bin_idx);
    if sum(valid_idx) == 0
         errordlg('Failed to discretize profile points.', 'Binning Error'); return;
    end
    
    % Fast calculation of mean values per bin (Forced to double to match NaN class)
    v_binned = accumarray(bin_idx(valid_idx), double(v_tube(valid_idx)), [bins, 1], @mean, NaN);
    if has_hgt
        h_binned = accumarray(bin_idx(valid_idx), double(h_tube(valid_idx)), [bins, 1], @mean, NaN);
    end

    % 8. Render Modern Dual-Axis Profile Plot
    h_prof_fig = figure;
    orient landscape;
    set(h_prof_fig, 'Name', ['Profile Plot ', momfig_name]);
    
    if has_hgt
        % --- Dual Axis: Velocity & Elevation ---
        yyaxis left
        % Plot raw scatter with alpha transparency (Modern Look)
        scatter(X_tube, v_tube, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        % Plot binned trend line
        plot(bin_centers, v_binned, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
        ylabel('Velocity (mm/yr)', 'Color', 'b', 'FontWeight', 'bold');
        set(gca, 'YColor', 'b');
        
        yyaxis right
        scatter(X_tube, h_tube, 15, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        plot(bin_centers, h_binned, 'k-s', 'LineWidth', 2, 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3]);
        ylabel('Elevation (m)', 'Color', [0.3 0.3 0.3], 'FontWeight', 'bold');
        set(gca, 'YColor', [0.3 0.3 0.3]);
        
    else
        % --- Single Axis: Velocity Only ---
        scatter(X_tube, v_tube, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        plot(bin_centers, v_binned, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
        ylabel('Velocity (mm/yr)', 'Color', 'b', 'FontWeight', 'bold');
    end
    
    hold off; grid on; set(gca, 'FontSize', 12);
    xlabel('Distance along profile (km)', 'FontWeight', 'bold');
    title(sprintf('Swath Profile (%d PS points, Radius: %d m)\nSource: %s', ...
                  n_pts_used, radiusfactor, momfig_name), 'Interpreter', 'none');
end