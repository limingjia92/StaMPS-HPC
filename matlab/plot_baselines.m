function [] = plot_baselines()
%PLOT_BASELINES Plot spatiotemporal baseline network for PS or SBAS
%   Unified function to automatically plot and export high-quality baseline 
%   networks for both Single-Master (PS) and Multi-Master (SBAS) modes.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          March 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   FEATURES:
%   - Auto-mode detection & direct topological indexing (No text file I/O).
%   - LSQ-based baseline reconstruction to prevent SBAS vector mismatches.
%   - Dynamic 8-digit date labels & auto-export to high-res PNG.
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, March 2011
%   ======================================================================

    disp('Generating Spatiotemporal Baseline Network Plot...');

    % 1. Load Data
    try
        load psver
        psname = ['ps', num2str(psver)];
        ps = load(psname);
    catch
        error('StaMPS-HPC: FileNotFound. Cannot find psver.mat or ps data. Please run step 1 first.');
    end

    % 2. Determine Processing Mode
    small_baseline_flag = getparm('small_baseline_flag');
    if isempty(small_baseline_flag)
        small_baseline_flag = 'n'; % Default to PS if not set
    end

    % 3. Setup Figure
    h_fig = figure('Name', 'Spatiotemporal Baseline Network', ...
                   'Position', [200, 200, 1000, 600], ... % Widened slightly for text labels
                   'Color', 'w');
    h_fig.WindowState = 'maximized';
    hold on;
    grid on;
    box on;

    % =====================================================================
    % MODE A: SBAS (Multi-Master)
    % =====================================================================
    if strcmpi(small_baseline_flag, 'y')
        disp('Mode: SBAS (Multi-Master Network)');
        
        if ~isfield(ps, 'ifgday_ix')
            error('StaMPS-HPC: Missing ''ifgday_ix'' in ps data. Cannot plot SBAS network.');
        end
        
        ifg_pairs = ps.ifgday_ix;
        n_ifg_total = size(ifg_pairs, 1);
        
        % Reconstruct Image Baselines from Interferogram Baselines
        if length(ps.bperp) == length(ps.day)
            image_bperp = double(ps.bperp);
        else
            G = zeros(n_ifg_total, length(ps.day));
            for i = 1:n_ifg_total
                G(i, ifg_pairs(i,1)) = -1; % Master image of the IFG
                G(i, ifg_pairs(i,2)) =  1; % Slave image of the IFG
            end
            
            % Fix Super Master to 0 to make the system solvable
            if isfield(ps, 'master_ix') && ~isempty(ps.master_ix)
                m_ix = ps.master_ix;
            else
                m_ix = 1; 
            end
            
            G_solve = G;
            G_solve(:, m_ix) = []; % Remove master column
            img_b_solve = G_solve \ double(ps.bperp); % Solve G * B_img = B_ifg
            
            image_bperp = zeros(length(ps.day), 1);
            idx_solve = setdiff(1:length(ps.day), m_ix);
            image_bperp(idx_solve) = img_b_solve;
        end
        
        % Filter out dropped interferograms (if any)
        drop_ifg_index = getparm('drop_ifg_index');
        if ~isempty(drop_ifg_index)
            keep_ix = setdiff(1:n_ifg_total, drop_ifg_index);
            ifg_pairs = ifg_pairs(keep_ix, :);
            disp(['   Note: Excluded ', num2str(length(drop_ifg_index)), ' dropped interferograms.']);
        end
        
        % Plot Network Edges (Lines)
        for i = 1:size(ifg_pairs, 1)
            idx1 = ifg_pairs(i, 1);
            idx2 = ifg_pairs(i, 2);
            plot([ps.day(idx1), ps.day(idx2)], [image_bperp(idx1), image_bperp(idx2)], ...
                 '-', 'Color', [0.3 0.6 0.8], 'LineWidth', 1.2);
        end
        
        if isfield(ps, 'master_ix') && ~isempty(ps.master_ix)
            m = ps.master_ix;
        else
            m = 1;
        end
        node_color = [1 0.8 0]; % Yellow for SBAS nodes

    % =====================================================================
    % MODE B: PS (Single-Master)
    % =====================================================================
    else
        disp('Mode: PS InSAR (Single-Master Network)');
        m = ps.master_ix;
        image_bperp = double(ps.bperp);
        
        % Plot Network Edges (Lines to Master)
        for i = 1:length(ps.day)
            if i ~= m
                plot([ps.day(m), ps.day(i)], [image_bperp(m), image_bperp(i)], ...
                     '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
            end
        end
        
        node_color = [0.3 0.6 0.8]; % Blue for PS nodes
    end

    % =====================================================================
    % 3.5 Unified Node Plotting & Date Labeling
    % =====================================================================
    % Plot Regular Network Nodes
    plot(ps.day, image_bperp, 'ko', 'MarkerFaceColor', node_color, 'MarkerSize', 7, 'LineWidth', 1);

    % Highlight Master Image
    plot(ps.day(m), image_bperp(m), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 14, 'LineWidth', 1.5);
    
    % Generate 8-digit date strings (e.g., '20220914')
    date_strs = cellstr(datestr(ps.day, 'yyyymmdd'));
    
    % Dynamic Y-offset to prevent text from covering the node
    y_offset = (max(image_bperp) - min(image_bperp)) * 0.02;
    if y_offset == 0, y_offset = 10; end
    
    % Add Date Labels to all nodes
    for i = 1:length(ps.day)
        if i == m
            % Master Date Label (Red, Bold, Slightly Larger)
            text(double(ps.day(i)), double(image_bperp(i) + y_offset), date_strs{i}, ...
                 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 14, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        else
            % Regular Date Label (Black)
            text(double(ps.day(i)), double(image_bperp(i) + y_offset), date_strs{i}, ...
                 'Color', 'k', 'FontSize', 12, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end

    % =====================================================================
    % 4. Aesthetics and Formatting
    % =====================================================================
    % Format X-axis as Dates
    x_ticks = get(gca, 'XTick');
    if length(x_ticks) >= 2 && (x_ticks(2) - x_ticks(1) < 60)
        date_fmt = 'yyyy-mm-dd';
    else
        date_fmt = 'yyyy-mm';
    end
    datetick('x', date_fmt, 'keepticks', 'keeplimits');  

    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    xlabel('Acquisition Date', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Perpendicular Baseline (m)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Spatiotemporal Baseline Network', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Ensure plot limits have some margins (Expand Top and Right for 45-deg text)
    y_lims = ylim;
    y_range = max(abs(y_lims));
    if y_range == 0, y_range = 100; end
    ylim([y_lims(1) - y_range*0.1 - 10, y_lims(2) + y_range*0.2 + 20]);
    
    x_lims = xlim;
    x_range = x_lims(2) - x_lims(1);
    xlim([x_lims(1) - x_range*0.05, x_lims(2) + x_range*0.1]);
    
    hold off;

    % =====================================================================
    % 5. Auto-Export to PNG
    % =====================================================================
    out_file = 'baseline.png';
    disp(['Saving high-resolution plot to: ', out_file]);
    
    % Use print for high-resolution 300 dpi output
    print(h_fig, '-dpng', '-r300', out_file);
    
end