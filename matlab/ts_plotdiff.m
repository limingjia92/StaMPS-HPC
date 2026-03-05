function ts_plotdiff(h_fig)
%TS_PLOTDIFF Interactive time series double difference plotting (HPC Optimized)
%   Triggered by the 'TS double diff.' UI button dynamically generated in ps_plot.
%
%   INPUTS:
%   h_fig    : Handle of the main scatter plot figure containing bound TS data.
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
%   1. Functionization & Data Binding: Encapsulated as a pure function using 
%      appdata, eliminating workspace pollution and fragile text-file passing.
%   2. Interactive Marker System: Implemented global tracking to seamlessly 
%      manage and clear multiple selection crosshairs across different tools.
%   3. SBAS Network Interception: Added strict dimensional safeguards to block 
%      invalid SM extractions on unlinked SBAS pairs, preventing crashes.
%   4. Visual Master Integration: Automatically injects the zero-deformation 
%      master date into both plots to ensure continuous and intuitive linear fits.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    % 1. Retrieve bound TS data file
    ts_matname = getappdata(h_fig, 'ts_matname');
    if isempty(ts_matname) || ~exist(ts_matname, 'file')
        errordlg('Cannot find bound TS data file. Please regenerate.', 'Data Error');
        return;
    end
    
    % 2. Get user-defined extraction radius
    mEditBox = findobj(h_fig, 'Tag', 'RadiusEditBox');
    if isempty(mEditBox)
        radiusfactor = 100;
    else
        radius_vals = str2num(char(get(mEditBox,'String'))); %#ok<ST2NM>
        if isempty(radius_vals), radiusfactor = 100; else, radiusfactor = radius_vals(1); end
    end
    
    % 3. Load exclusive TS data for this figure
    ts_data = load(ts_matname);
    ph_mm = ts_data.ph_mm;
    lonlat = ts_data.lonlat;
    day = ts_data.day;
    master_day = ts_data.master_day;
    momfig_name = ['(', get(h_fig, 'Name'), ')'];
    
    % 4. Activate main figure and await dual user selections
    figure(h_fig);
    
    % Clear previous residual markers
    old_marker = getappdata(h_fig, 'ts_marker_handle');
    if ~isempty(old_marker) && all(isgraphics(old_marker))
        delete(old_marker); 
    end

    hold(gca, 'on');
    disp(['Select 1ST point (Reference) for TS Double Diff (Radius: ', num2str(radiusfactor), ' m)...']);
    try
        [lon1, lat1] = ginput(1);
        h_marker1 = plot(gca, lon1, lat1, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        
        % Save 1st marker immediately to prevent ghosting if 2nd selection is cancelled
        setappdata(h_fig, 'ts_marker_handle', h_marker1); 
        disp(['  -> Pt 1 Coordinates: ', num2str(lon1), ', ', num2str(lat1)]);
    catch
        disp('Point 1 selection cancelled.'); hold(gca, 'off'); return;
    end
    
    disp(['Select 2ND point (Target) for TS Double Diff (Radius: ', num2str(radiusfactor), ' m)...']);
    try
        [lon2, lat2] = ginput(1);
        h_marker2 = plot(gca, lon2, lat2, 'bx', 'MarkerSize', 12, 'LineWidth', 2.5); 
        
        % Update appdata with both marker handles
        setappdata(h_fig, 'ts_marker_handle', [h_marker1, h_marker2]); 
        disp(['  -> Pt 2 Coordinates: ', num2str(lon2), ', ', num2str(lat2)]);
    catch
        disp('Point 2 selection cancelled.'); hold(gca, 'off'); return;
    end
    hold(gca, 'off');

    % 5. Find valid PS points within search radii
    search_radius_deg = radiusfactor / 1000 / 111;
    
    % Locate Point 1 (Reference)
    dist1 = distance(lonlat(:,2), lonlat(:,1), lat1, lon1);
    pts_near_ix1 = dist1 <= search_radius_deg;
    n_pts_near1 = sum(pts_near_ix1);
    if n_pts_near1 == 0
        [~, min_idx1] = min(dist1); pts_near_ix1(min_idx1) = true; n_pts_near1 = 1;
        disp('Warning: No points within Pt 1 radius. Selecting the closest point.');
    end
    
    % Locate Point 2 (Target)
    dist2 = distance(lonlat(:,2), lonlat(:,1), lat2, lon2);
    pts_near_ix2 = dist2 <= search_radius_deg;
    n_pts_near2 = sum(pts_near_ix2);
    if n_pts_near2 == 0
        [~, min_idx2] = min(dist2); pts_near_ix2(min_idx2) = true; n_pts_near2 = 1;
        disp('Warning: No points within Pt 2 radius. Selecting the closest point.');
    end
    
    % 6. Extract and compute double difference (Pt 1 - Pt 2)
    ts1 = mean(ph_mm(pts_near_ix1, :), 1)';
    ts2 = mean(ph_mm(pts_near_ix2, :), 1)';
    ts_diff = ts1 - ts2;
    
    % 7. SBAS invalid network interception
    if length(ts_diff) ~= length(day)
        err_msg = { 'Time series dimension mismatch (Physical logic conflict)!', '', ...
                    'Cause: You are using lowercase ''v'' in a Small Baseline (SBAS) network. This prevents time-series differencing.', ...
                    'Solution: Close this window and use uppercase ''V'' to force SBAS-to-Single-Master inversion.' };
        errordlg(err_msg, 'StaMPS-HPC Safety Guard');
        return;
    end
    
    % 8. Linear rate fitting (relative velocity)
    G = [ones(size(day)), day - master_day]; 
    m = lscov(G, ts_diff); 
    
    % 9. Visual enhancement: Auto-append master zero point
    if ~ismember(master_day, day)
        day_plot = [day; master_day];
        ts_diff_plot = [ts_diff; 0]; % Diff at master is exactly 0
        [day_plot, sort_idx] = sort(day_plot);
        ts_diff_plot = ts_diff_plot(sort_idx);
    else
        day_plot = day;
        ts_diff_plot = ts_diff;
    end
    
    % Recalculate fitted line including zero point for plotting
    G_plot = [ones(size(day_plot)), day_plot - master_day];
    ts_hat_diff = G_plot * m;
    
    % 10. Render location context map
    h_loc_fig = figure;
    set(h_loc_fig, 'Name', 'Selection Context Map');
    
    plot(lonlat(pts_near_ix1, 1), lonlat(pts_near_ix1, 2), 'sr', 'MarkerFaceColor', 'r'); hold on;
    plot(lonlat(pts_near_ix2, 1), lonlat(pts_near_ix2, 2), 'sb', 'MarkerFaceColor', 'b');
    plot(lon1, lat1, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 15);
    plot(lon2, lat2, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 15);
    legend('Pt 1 (Reference)', 'Pt 2 (Target)', 'Click Center 1', 'Click Center 2', 'Location', 'Best');
    title(sprintf('Location Context (Radius: %d m)', radiusfactor));
    axis image; grid on; hold off;
    
    % 11. Render modernized TS Double Diff plot
    h_ts_fig = figure;
    orient landscape;
    set(h_ts_fig, 'Name', ['TS Double Diff ', momfig_name]);
    
    DAY = datetime(day_plot, 'ConvertFrom', 'datenum');
    
    % Plot fitted line and scatter points
    plot(DAY, ts_hat_diff, '-r', 'LineWidth', 1.5); hold on;
    plot(DAY, ts_diff_plot, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    
    % Highlight master date
    master_dt = datetime(master_day, 'ConvertFrom', 'datenum');
    plot(master_dt, 0, 'pg', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    
    hold off; grid on; set(gca, 'FontSize', 12);
    ylabel('LOS (mm)', 'FontWeight', 'bold');
    title(sprintf('TS Double Difference\nPt 1 (%d pts) - Pt 2 (%d pts)\nSource: %s', ...
                  n_pts_near1, n_pts_near2, momfig_name), 'Interpreter', 'none');
end