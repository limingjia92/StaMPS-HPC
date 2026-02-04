function [] = sb_load_initial_gamma(endian)
% SB_LOAD_INITIAL_GAMMA (HPC Optimized Version)
%   Initializes the SBAS workflow by loading files using GAMMA outputs.
%
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Refactored by: Mingjia Li
%   Date:          January 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization & Logic Updates:
%   1. Dependency Consolidation: Merged 'readparm' into internal function.
%   2. Logic Update: Optimized Heading/Incidence angle reading priority
%      (File -> Calculation) and ensured correct save filenames (la1, inc1, head1).
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, Feb 2013
%   ======================================================================

    if nargin < 1
       endian = 'b'; % Default to big-endian for GAMMA
    end

    logit;
    logit('Loading data into matlab from GAMMA (SBAS)...\n');

    %% --- 0. Configuration & File Paths ---
    phname   = './pscands.1.ph';    % Complex phase
    ijname   = './pscands.1.ij';    % ID, Azimuth, Range
    llname   = './pscands.1.ll';    % Lon, Lat
    hgtname  = './pscands.1.hgt';   % Height
    daname   = './pscands.1.da';    % Amplitude Dispersion
    headname = './pscands.1.head';  % Heading angle file
    incname  = './pscands.1.inc';   % Incidence angle file
    rscname  = '../rsc.txt';        % Resource file pointing to SLC par
    pscname  = '../pscphase.in';    % Interferogram list

    psver = 1;

    %% --- 1. Read RSLC Parameters ---
    if ~exist(rscname, 'file')
        error([rscname, ' does not exist']);
    end
    fid = fopen(rscname, 'r');
    rslcpar_cell = textscan(fid, '%s');
    fclose(fid);
    rslcpar = rslcpar_cell{1}{1}; 

    if ~exist(rslcpar, 'file')
        error(['RSLC parameter file not found: ', rslcpar]);
    end

    master_day_val = str2double(rslcpar(end-16:end-9));
    year = floor(master_day_val / 10000);
    month = floor((master_day_val - year * 10000) / 100);
    day_val = master_day_val - year * 10000 - month * 100;
    master_day = datenum(year, month, day_val);

    %% --- 2. Read Interferogram List & Dates (SBAS Logic) ---
    if ~exist(pscname, 'file')
        error([pscname, ' does not exist']);
    end
    fid = fopen(pscname, 'r');
    ifgs_raw = textscan(fid, '%s');
    fclose(fid);
    
    ifgs = ifgs_raw{1}(2:end); % Skip header
    nb = length(ifgs{1});      % Length of filename string
    n_ifg = length(ifgs);
    
    ifgday = zeros(n_ifg, 2);
    for i = 1:n_ifg
        ifgday(i, 1) = str2double(ifgs{i}(nb-21:nb-14)); % Primary
        ifgday(i, 2) = str2double(ifgs{i}(nb-12:nb-5));  % Secondary
    end

    % Convert to datenums
    year = floor(ifgday / 10000);
    month = floor((ifgday - year * 10000) / 100);
    day_val = ifgday - year * 10000 - month * 100;
    ifgday_dn = datenum(year, month, day_val);
    
    % Unique days list
    [day, ~, ifgday_ix] = unique(ifgday_dn(:));
    ifgday_ix = reshape(ifgday_ix, n_ifg, 2);

    ifgday = ifgday_dn;

    n_image = length(day);
    master_ix = sum(day < master_day) + 1;

    %% --- 3. Radar Parameters & Platform ---
    freq = get_gamma_par(rslcpar, 'radar_frequency:', 1);
    lambda = 299792458 / freq;
    setparm('lambda', lambda, 1);

    sensor = get_gamma_par(rslcpar, 'sensor:', 1);
    if contains(sensor, 'ASAR')
        platform = 'ENVISAT';
    else
        platform=sensor; % S1A for Sentinel-1A
    end
    setparm('platform', platform, 1);

    %% --- 4. Geometry: Heading & Incidence (Logic Updated) ---
    % -- Heading --
    if exist(headname, 'file')
        logit(['Reading heading from ' headname]);
        fid = fopen(headname, 'r');
        head = fread(fid, [1, inf], 'float', endian)';
        fclose(fid);
        heading_angle = mean(head);
        head = head ./180 * pi; % convert from degree to radius
    else
        logit('Calculating heading from RSLC parameters...');
        heading_angle = get_gamma_par(rslcpar, 'heading:', 1);
        head = []; 
    end
    setparm('heading', heading_angle, 1);

    % -- Geometry Calculations --
    ij = load(ijname);
    n_ps = size(ij, 1);
    
    rps = get_gamma_par(rslcpar, 'range_pixel_spacing', 1);
    rgn = get_gamma_par(rslcpar, 'near_range_slc', 1);
    se  = get_gamma_par(rslcpar, 'sar_to_earth_center', 1);
    re  = get_gamma_par(rslcpar, 'earth_radius_below_sensor', 1);
    rgc = get_gamma_par(rslcpar, 'center_range_slc', 1);

    rg = rgn + ij(:,3) * rps; % Range distance
    mean_range = rgc; % mean range distance
    
    % -- Look Angle (Calculated) --
    % Always calculated to ensure consistency, stored in la1
    la = acos((se^2 + rg.^2 - re^2) ./ (2 * se * rg));

    % -- Incidence Angle --
    if exist(incname, 'file')
        logit(['Reading incidence from ' incname]);
        fid = fopen(incname, 'r');
        inc = fread(fid, [1, inf], 'float', endian)';
        fclose(fid);
        inc = inc ./180 * pi; % convert from degree to radius
        mean_incidence = mean(inc);
    else
        logit('Calculating incidence from geometry...');
        inc = acos((se^2 - re^2 - rg.^2) ./ (2 * re * rg));
        mean_incidence = mean(inc);
    end

    %% --- 5. Baseline Calculation ---
    naz = get_gamma_par(rslcpar, 'azimuth_lines', 1);
    prf = get_gamma_par(rslcpar, 'prf', 1);
    mean_az = naz / 2 - 0.5;

    bperp_mat = zeros(n_ps, n_ifg, 'single');
    
    logit('Calculating perpendicular baselines...');
    for i = 1:n_ifg
        basename = [ifgs{i}(1:nb-4), 'base'];
        if exist(basename, 'file')
            B_TCN = get_gamma_par(basename, 'initial_baseline(TCN):', 3, 0);
            BR_TCN = get_gamma_par(basename, 'initial_baseline_rate:', 3, 0);
            
            if iscell(B_TCN), B_TCN = [B_TCN{:}]; end
            if iscell(BR_TCN), BR_TCN = [BR_TCN{:}]; end
            
            bc = B_TCN(2) + BR_TCN(2) * (ij(:,2) - mean_az) / prf;
            bn = B_TCN(3) + BR_TCN(3) * (ij(:,2) - mean_az) / prf;
            
            % Convert (T)CN to perp-para
            bperp_mat(:, i) = bc .* cos(la) - bn .* sin(la);
        else
            fprintf('Warning: Baseline file %s not found. Setting to 0.\n', basename);
        end
    end
    bperp = mean(bperp_mat)';

    %% --- 6. Load Phase Data (Vectorized) ---
    logit('Reading phase data...');
    fid = fopen(phname, 'r', endian);
    % Vectorized read: [2*n_ps, n_ifg] -> Faster than loop
    raw_ph = fread(fid, [2 * n_ps, n_ifg], 'float=>single');
    fclose(fid);
    
    % Convert to complex [n_ps, n_ifg]
    ph = complex(raw_ph(1:2:end, :), raw_ph(2:2:end, :));
    clear raw_ph;

    %% --- 7. Load Lon/Lat & Coordinate Rotation ---
    if ~exist(llname, 'file')
        error([llname, ' does not exist']);
    end
    fid = fopen(llname, 'r');
    lonlat = fread(fid, [2, inf], 'float', endian)';
    fclose(fid);

    ll0 = (max(lonlat) + min(lonlat)) / 2;
    xy = llh2local(lonlat', ll0)' * 1000;

    % Rotation
    theta = (180 - heading_angle) * pi / 180;
    if theta > pi, theta = theta - 2 * pi; end
    
    rotm = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    xynew = (rotm * xy')'; 
    
    range_old = max(xy) - min(xy);
    range_new = max(xynew) - min(xynew);
    
    if range_new(1) < range_old(1) && range_new(2) < range_old(2)
        xy = xynew;
        logit(['Rotating by ', num2str(theta * 180 / pi), ' degrees']);
    end
    
    xy = single(xy);
    [~, sort_ix] = sortrows(xy, [2, 1]); 
    xy = xy(sort_ix, :);
    
    xy = [[1:n_ps]', xy];
    xy(:, 2:3) = round(xy(:, 2:3) * 1000) / 1000; 

    % Apply Sort
    ph = ph(sort_ix, :);
    ij = ij(sort_ix, :);
    ij(:, 1) = 1:n_ps; 
    lonlat = lonlat(sort_ix, :);
    bperp_mat = bperp_mat(sort_ix, :);
    la = la(sort_ix);
    inc = inc(sort_ix);
    if ~isempty(head), head = head(sort_ix); end

    %% --- 8. NaN Filtering & Saving ---
    % Robust NaN filtering (Consistent with PS optimized version)
    ix_nan = (sum(isnan(lonlat), 2) >= 1 | sum(isnan(ph), 2) >= 1);
    
    ij(ix_nan, :) = [];
    lonlat(ix_nan, :) = [];
    xy(ix_nan, :) = [];
    ph(ix_nan, :) = [];
    bperp_mat(ix_nan, :) = [];
    la(ix_nan, :) = [];
    inc(ix_nan, :) = [];
    if ~isempty(head), head(ix_nan, :) = []; end

    n_ps = size(lonlat, 1);
    ij(:, 1) = [1:n_ps]';
    xy(:, 1) = [1:n_ps]';

    % -- Main SBAS Save --
    savename = ['ps', num2str(psver)];
    stamps_save(savename, ij, lonlat, xy, bperp, day, master_day, master_ix, ifgday, ifgday_ix, n_ifg, n_image, n_ps, sort_ix, ll0, master_ix, mean_incidence, mean_range);
    save('psver', 'psver');

    % -- Phase Save --
    phsavename = ['ph', num2str(psver)];
    stamps_save(phsavename, ph);

    % -- Baseline Grid Save --
    bpsavename = ['bp', num2str(psver)];
    stamps_save(bpsavename, bperp_mat);

    % -- Look Angle Save (la1) --
    lasavename = ['la', num2str(psver)];
    stamps_save(lasavename, la);

    % -- Incidence Angle Save (inc1) --
    incsavename = ['inc', num2str(psver)];
    stamps_save(incsavename, inc);

    % -- Heading Save (head1) --
    if ~isempty(head)
        headsavename = ['head', num2str(psver)];
        stamps_save(headsavename, head);
    end

    %% --- 9. Auxiliary Files (DA / HGT) ---
    if exist(daname, 'file')
        D_A = load(daname);
        D_A = D_A(sort_ix);
        D_A(ix_nan, :) = [];
        stamps_save(['da', num2str(psver)], D_A);
    end

    if exist(hgtname, 'file')
        fid_hgt = fopen(hgtname, 'r');
        hgt = fread(fid_hgt, [1, inf], 'float', endian)';
        fclose(fid_hgt);
        hgt = hgt(sort_ix);
        hgt(ix_nan, :) = [];
        stamps_save(['hgt', num2str(psver)], hgt);
    end

    logit(1);
end

%% --- Internal Helper Functions ---
function [values] = get_gamma_par(fname, parm, numval, log_flag)
% GET_GAMMA_PAR Reads a parameter from a GAMMA file
% Replacement for the external 'readparm.m' 
% 
% Inputs:
%   fname: File path
%   parm:  Parameter string key (e.g., 'radar_frequency:')
%   numval: Number of values to read after the key
%   log_flag
%   log_flag: Boolean to print log

    if nargin < 3, numval = 1; end
    if nargin < 4, log_flag = 1; end

    txt = fileread(fname);
    idx = strfind(txt, parm);
    if isempty(idx)
        error(['Parameter ' parm ' not found in ' fname]);
    end
    
    curr_pos = idx(1) + length(parm);
    txt_len = length(txt);
    
    while curr_pos <= txt_len
        char_now = txt(curr_pos);
        if isspace(char_now) || strcmp(char_now, ':')
            curr_pos = curr_pos + 1;
        else
            break; 
        end
    end
    
    sub_txt = txt(curr_pos:end);
    data = textscan(sub_txt, '%s', numval);
    tokens = data{1};
    
    if isempty(tokens)
        values = [];
        return;
    end

    val_cell = cell(1, numval);
    log_str = [parm, ' '];
    
    for i = 1:numval
        token = tokens{i};
        num = str2double(token);
        
        if isnan(num)
            val_cell{i} = token; 
        else
            val_cell{i} = num;   
        end
        log_str = [log_str, token, ' '];
    end
    
    if numval == 1
        values = val_cell{1}; 
    else
        % If all elements are numbers, return a vector instead of cell
        if all(cellfun(@isnumeric, val_cell))
            values = [val_cell{:}];
        else
            values = val_cell;
        end
    end

    if log_flag
        logit(log_str);
    end
end