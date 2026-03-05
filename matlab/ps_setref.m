function [ref_ps] = ps_setref(ps2)
%PS_SETREF Select reference pixels based on location (HPC Optimized)
%   Selects reference PS candidates based on geographic bounding box or 
%   circular radius from a center point.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Framework Optimization: Enabled direct structure input to minimize I/O overhead.
%   2. Algorithm Optimization: Vectorized distance calculations for radius filtering.
%   3. Logic Cleanup: Removed deprecated 'ref_x/ref_y' legacy filtering logic.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    % =====================================================================
    % 1. Input Handling & Data Extraction
    % =====================================================================
    if nargin < 1
        % Legacy Mode: Load data from disk if not provided
        load psver
        psname = ['ps', num2str(psver)];
        fprintf('PS_SETREF: No inputs provided. Loading %s...\n', psname);
        ps2 = load(psname);
    end
    
    % Ensure required fields exist
    if ~isfield(ps2, 'lonlat')
        error('PS_SETREF: Input structure must contain ''lonlat'' field.');
    end
    
    % Extract master reference origin (ll0)
    if isfield(ps2, 'll0')
        ll0 = ps2.ll0;
    else
        error('PS_SETREF: Cannot find ll0 parameter.');
    end

    lonlat = ps2.lonlat;
    n_ps = size(lonlat, 1);

    % =====================================================================
    % 2. Parameter Retrieval
    % =====================================================================
    ref_lon = getparm('ref_lon');
    ref_lat = getparm('ref_lat');
    ref_centre_lonlat = getparm('ref_centre_lonlat');
    ref_radius = getparm('ref_radius');

    % =====================================================================
    % 3. Filtering Logic
    % =====================================================================
    
    % --- Filter A: Geographic Bounding Box ---
    mask_box = true(n_ps, 1);
    
    if ~isempty(ref_lon) && ~isempty(ref_lat)
        mask_box = lonlat(:,1) > ref_lon(1) & lonlat(:,1) < ref_lon(2) & ...
                   lonlat(:,2) > ref_lat(1) & lonlat(:,2) < ref_lat(2);
    end

    % --- Filter B: Circular Radius Search ---
    mask_radius = true(n_ps, 1);
    
    if ~isempty(ref_radius) && ref_radius < inf && ref_radius > 0
        % Convert Reference Center to local coordinates (m)
        ref_xy_center = llh2local(ref_centre_lonlat', ll0) * 1000;

        % Convert Candidates to local coordinates (m) relative to the SAME origin (ll0)
        % Note: This ensures geometric consistency regardless of original ps2.xy 
        xy_local = llh2local(lonlat', ll0)' * 1000;

        % Vectorized Squared Euclidean Distance Calculation
        dist_sq = (xy_local(:,1) - ref_xy_center(1)).^2 + (xy_local(:,2) - ref_xy_center(2)).^2;
        
        mask_radius = dist_sq <= ref_radius^2;
        
    elseif ref_radius == -inf
        % Explicitly disable reference selection
        mask_radius = false(n_ps, 1);
    end

    % =====================================================================
    % 4. Selection & Reporting
    % =====================================================================
    
    % Intersect filters
    ref_ps = find(mask_box & mask_radius);
    
    % Fallback: If no points selected (and not explicitly disabled), select ALL
    if isempty(ref_ps) && ~(exist('ref_radius','var') && ref_radius == -inf)
        fprintf('PS_SETREF: No specific points selected. Defaulting to ALL points.\n');
        ref_ps = (1:n_ps)';
    else
        fprintf('PS_SETREF: Selected %d reference points.\n', length(ref_ps));
    end

end