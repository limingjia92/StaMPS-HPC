function ts_plot(h_fig)
%TS_PLOT Interactive time series extraction and plotting (HPC Optimized)
%   Triggered by the UI button dynamically generated in ps_plot.
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
%   1. Interactive UX Enhancements: Implemented global marker tracking via appdata 
%      to dynamically clear and render selection crosshairs without base map redraws.
%   2. Robust Data Binding: Replaced slow global variables with direct figure-bound 
%      appdata for encapsulated, conflict-free memory access.
%   3. SBAS Network Interception: Added a strict dimensional safeguard to prevent 
%      mathematical collapse when attempting SM extraction on unlinked SBAS pairs.
%   4. Visual Master Integration: Automatically injects the zero-deformation master 
%      date into the plotting arrays to ensure a continuous and intuitive linear fit.
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
    
    % 4. Activate main figure and await user selection
    figure(h_fig);
    
    % Clear previous residual marker
    old_marker = getappdata(h_fig, 'ts_marker_handle');
    if ~isempty(old_marker) && all(isgraphics(old_marker))
        delete(old_marker); 
    end
    
    disp(['Select point for Time Series extraction (Radius: ', num2str(radiusfactor), ' m)...']);
    try
        % Capture mouse click on the axes
        [lon_req, lat_req] = ginput(1);
        disp(['  -> Point Coordinates: ', num2str(lon_req), ', ', num2str(lat_req)]);
        
        % Draw a prominent red cross marker
        hold(gca, 'on');
        h_marker = plot(gca, lon_req, lat_req, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        
        % Save marker handle for future deletion
        setappdata(h_fig, 'ts_marker_handle', h_marker);
        hold(gca, 'off');
    catch
        disp('Point selection cancelled.');
        return;
    end
    
    % 5. Find valid PS points within search radius
    dist = distance(lonlat(:,2), lonlat(:,1), lat_req, lon_req);
    search_radius_deg = radiusfactor / 1000 / 111;
    pts_near_ix = dist <= search_radius_deg;
    n_pts_near = sum(pts_near_ix);
    
    if n_pts_near == 0
        [~, min_index] = min(dist);
        pts_near_ix(min_index) = true;
        n_pts_near = 1;
        disp('Warning: No points within radius. Selecting the closest single point instead.');
    end
    
    % 6. Extract deformation data
    ts_matrix = ph_mm(pts_near_ix, :);
    ts = mean(ts_matrix, 1)';
    
    % 7. SBAS invalid network interception
    if length(ts) ~= length(day)
        err_msg = { 'Time series dimension mismatch (Physical logic conflict)!', '', ...
                    sprintf('Point has %d phase values, but time axis has %d dates.', length(ts), length(day)), ...
                    'Cause: You are using lowercase ''v'' in a Small Baseline (SBAS) network. This extracts raw interferogram pairs, which cannot be plotted on a single timeline.', ...
                    'Solution: Close this window and use uppercase ''V'' (e.g., ''V-da(era5)o'') to force SBAS-to-Single-Master inversion.' };
        errordlg(err_msg, 'StaMPS-HPC Safety Guard');
        return;
    end
    
    % 8. Linear rate fitting (excluding master 0 point for mathematical rigor)
    G = [ones(size(day)), day - master_day]; 
    m = lscov(G, ts); 
    
    % 9. Visual enhancement: Auto-append master zero point
    if ~ismember(master_day, day)
        day_plot = [day; master_day];
        ts_plot_vals = [ts; 0];
        [day_plot, sort_idx] = sort(day_plot);
        ts_plot_vals = ts_plot_vals(sort_idx);
    else
        day_plot = day;
        ts_plot_vals = ts;
    end
    
    % Recalculate fitted line including zero point for plotting
    G_plot = [ones(size(day_plot)), day_plot - master_day];
    ts_hat = G_plot * m;
    
    % 10. Render modernized TS plot
    h_ts_fig = figure;
    orient landscape;
    set(h_ts_fig, 'Name', ['Time series plot for ', num2str(n_pts_near), ' point(s) ', momfig_name]);
    
    DAY = datetime(day_plot, 'ConvertFrom', 'datenum');
    
    % Plot fitted line and scatter points
    plot(DAY, ts_hat, '-r', 'LineWidth', 1.5); hold on;
    plot(DAY, ts_plot_vals, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    
    % Highlight master date
    master_dt = datetime(master_day, 'ConvertFrom', 'datenum');
    plot(master_dt, 0, 'pg', 'MarkerFaceColor', 'g', 'MarkerSize', 10); 
    
    hold off; grid on; set(gca, 'FontSize', 12);
    ylabel('LOS (mm)', 'FontWeight', 'bold');
    title(sprintf('TS Plot: %d PS points within %dm radius\nSource: %s', ...
                  n_pts_near, radiusfactor, momfig_name), 'Interpreter', 'none');
end