function [] = ps_load_initial_isce
% PS_LOAD_INITIAL_ISCE (HPC Optimized Version)
%   Initializes the PS workflow by load files using ISCE outputs
%   
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Dependency Consolidation: Merged auxiliary functions to reduce file handle overhead.
%   2. Vectorized I/O: Replaced iterative loops with matrix `fread` and fast Regex parsing.
%   3. Algorithmic Efficiency: Substituted `griddata` with `interp2` for regular grid interpolation.
%   4. Add matfile head1.mat for save heading angle.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

    logit;
    logit('Loading data into matlab from ISCE (PS)...\n');

    % configuration parameter
    phname = 'pscands.1.ph';      % for each PS candidate, a float complex value for each ifg
    ijname = 'pscands.1.ij';      % ID# Azimuth# Range# 1 line per PS candidate
    bperpname = 'bperp.1.in';     % in meters 1 line per ifg
    dayname = 'day.1.in';         % YYYYMMDD, 1 line per ifg
    masterdayname = 'reference_day.1.in';  % YYYYMMDD
    llname = 'pscands.1.ll';      % 2 float values (lon and lat) per PS candidate
    daname = 'pscands.1.da';      % 1 float value per PS candidate
    hgtname = 'pscands.1.hgt';    % 1 float value per PS candidate
    incname = 'inc_angle.raw';    % grid of incidence angle values
    headname = 'heading.raw';     % satellite heading
    headingname = 'heading.1.in'; % mean satellite heading    
    lambdaname = 'lambda.1.in';   % wavelength
    calname = 'calamp.out';       % amplitide calibrations
    widthname = 'width.txt';      % width of interferograms
    lenname = 'len.txt';          % length of interferograms

    psver = 1;
    incsavename = ['inc', num2str(psver)];
    headsavename = ['head', num2str(psver)];

    % auxiliary function
    fix_path = @(x) conditional_path(x);

    %% --- 1. reading data information ---
    dayname = fix_path(dayname);
    day_data = load(dayname);
    year = floor(day_data / 10000);
    month = floor((day_data - year * 10000) / 100);
    monthday = day_data - year * 10000 - month * 100;
    slave_day = datenum(year, month, monthday);
    [slave_day, day_ix] = sort(slave_day);

    masterdayname = fix_path(masterdayname);
    master_day_val = load(masterdayname);
    master_day_yyyymmdd = master_day_val;
    myear = floor(master_day_val / 10000);
    mmonth = floor((master_day_val - myear * 10000) / 100);
    mday = master_day_val - myear * 10000 - mmonth * 100;
    master_day = datenum(myear, mmonth, mday);

    master_ix = sum(slave_day < master_day) + 1;
    day = [slave_day(1:master_ix-1); master_day; slave_day(master_ix:end)];

    %% --- 2. perpendicular baseline ---
    bperpname = fix_path(bperpname);
    bperp = load(bperpname);
    bperp = bperp(day_ix);
    bperp = [bperp(1:master_ix-1); 0; bperp(master_ix:end)];
    n_ifg = size(bperp, 1);
    n_image = n_ifg;

    %% --- 3. Heading & Lambda ---
    headingname = fix_path(headingname);
    heading = load(headingname);
    if isempty(heading), error('heading.1.in is empty'); end
    setparm('heading', heading, 1);

    lambdaname = fix_path(lambdaname);
    lambda = load(lambdaname);
    setparm('lambda', lambda, 1);
    
    if abs(lambda - 0.0554658) < 10e-4
        platform='S1A';
    elseif abs(lambda - 0.0562356) < 5e-4
        platform='ENVISAT';
    elseif abs(lambda - 0.236057) < 1e2
        platform = 'ALOS';
    elseif abs(lambda - 0.0310666) < 5e-4
        platform = 'TERRASAR';
    else
        platform = '';
    end
    setparm('platform', platform, 1);
    
    %% --- 4. Load PS Indices and Amplitude Calibration ---
    ij = load(ijname);
    n_ps = size(ij, 1);

    calname = fix_path(calname);
    if exist(calname, 'file')
        fid_cal = fopen(calname);
        cal_raw = textscan(fid_cal, '%s %f');
        fclose(fid_cal);
        calfile = cal_raw{1};
        calvals_raw = cal_raw{2};
        
        caldate = zeros(length(calfile), 1);
        for i = 1:length(calfile)
            date_pattern = regexp(calfile{i}, '\d{8}', 'match');
            if ~isempty(date_pattern)
                caldate(i) = str2double(date_pattern{end}); 
            else
                if contains(lower(calfile{i}), 'reference') || contains(lower(calfile{i}), 'master')
                    caldate(i) = master_day_yyyymmdd; % 
                else
                    caldate(i) = -1; %
                end
            end
        end
        
        slave_dates_sorted = day_data(day_ix); 
        calconst = ones(length(slave_dates_sorted), 1);
        for k = 1:length(slave_dates_sorted)
            target_date = slave_dates_sorted(k);
            
            idx = find(caldate == target_date);
            
            if ~isempty(idx)
                calconst(k) = calvals_raw(idx(1));
            else
                fprintf('Warning: Calibration constant for date %d not found in %s. Using 1.0.\n', target_date, calname);
            end
        end

    else
        calconst = ones(length(day_ix), 1);
    end

    %% --- 5. Phase data ---
    fid_ph = fopen(phname, 'r');
    if fid_ph == -1, error(['Cannot open ' phname]); end
    
    %  [2 * n_ps, n_ifg-1]
    raw_ph_data = fread(fid_ph, [n_ps * 2, n_ifg - 1], 'float=>single'); 
    fclose(fid_ph);

    ph = complex(raw_ph_data(1:2:end, :), raw_ph_data(2:2:end, :));
    clear raw_ph_data; % release memory

    ph = ph(:, day_ix);
    ph = bsxfun(@rdivide, ph, calconst'); % scale amplitudes
    % insect reference image phase as 1
    ph = [ph(:, 1:master_ix-1), ones(n_ps, 1, 'single'), ph(:, master_ix:end)];

    %% --- 6. Lon & LAT ---
    fid_ll = fopen(llname, 'r');
    lonlat = fread(fid_ll, [2, inf], 'float');
    lonlat = lonlat';
    fclose(fid_ll);

    ll0 = (max(lonlat) + min(lonlat)) / 2;
    xy = llh2local(lonlat', ll0) * 1000;
    xy = xy';

    % coordinate rotate
    theta = (180 - heading) * pi / 180;
    if theta > pi, theta = theta - 2 * pi; end
    rotm = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    xynew = (rotm * xy')'; 
    
    range_old = max(xy) - min(xy);
    range_new = max(xynew) - min(xynew);
    if range_new(1) < range_old(1) && range_new(2) < range_old(2)
        xy = xynew;
        fprintf('Rotating by %.2f degrees\n', theta * 180 / pi);
    end
    
    xy = single(xy);
    [~, sort_ix] = sortrows(xy, [2, 1]); % sort in ascending y order
    xy = xy(sort_ix, :);
    xy = [[1:n_ps]', xy];
    xy(:, 2:3) = round(xy(:, 2:3) * 1000) / 1000; % round to mm

    ph = ph(sort_ix, :);
    ij = ij(sort_ix, :);
    ij(:, 1) = 1:n_ps; % reset ID
    lonlat = lonlat(sort_ix, :);

    %% --- 7. save matfile ---
    savename = ['ps', num2str(psver)];
    stamps_save(savename, ij, lonlat, xy, bperp, day, master_day, master_ix, n_ifg, n_image, n_ps, sort_ix, ll0, calconst, master_ix, day_ix);
    save('psver', 'psver');

    phsavename = ['ph', num2str(psver)];
    stamps_save(phsavename, ph);

    %% --- 8. processing DA and HGT (if exist) ---
    if exist(daname, 'file')
        D_A = load(daname);
        D_A = D_A(sort_ix);
        stamps_save(['da', num2str(psver)], D_A);
    end

    if exist(hgtname, 'file')
        fid_hgt = fopen(hgtname, 'r');
        hgt = fread(fid_hgt, [1, inf], 'float')';
        fclose(fid_hgt);
        hgt = hgt(sort_ix);
        stamps_save(['hgt', num2str(psver)], hgt);
    end

    %% --- 9. read Width/Len ---
    widthname = fix_path(widthname);
    width = load(widthname);
    lenname = fix_path(lenname);
    len = load(lenname);

    %% --- 10. Processing Incidence/ Heading ---
    target_geo_file = '';
    target_save_name = '';

    incname = search_file_upwards(incname);
    
    if exist(incname, 'file')
        target_geo_file = incname;
        target_save_name = incsavename;
    end
    
    if ~isempty(target_geo_file)
        grid_data = load_isce_internal(target_geo_file);
    end

    if exist('grid_data', 'var')
        idx_sub = sub2ind(size(grid_data), ij(:,2) + 1, ij(:,3) + 1);
        inc = grid_data(idx_sub);
        inc = inc ./180 * pi; % convert from degree to radius
        stamps_save(target_save_name, inc);
        mean_incidence = mean(inc); % append mean_incidence value
        savename = ['ps', num2str(psver)];
        stamps_save(savename, mean_incidence);
        clear grid_data;
    end

    target_geo_file = '';
    target_save_name = '';

    headname = search_file_upwards(headname);
    
    if exist(headname, 'file')
        target_geo_file = headname;
        target_save_name = headsavename;
    end
    
    if ~isempty(target_geo_file)
        grid_data = load_isce_internal(target_geo_file);
    end

    if exist('grid_data', 'var')
        idx_sub = sub2ind(size(grid_data), ij(:,2) + 1, ij(:,3) + 1);
        head = grid_data(idx_sub);
        head = head * pi ./ 180; % convert from degree to radius
        stamps_save(target_save_name, head);
        clear grid_data;
    end

    %% --- 11. processing Baseline Grid  ---
    bpsavename = ['bp', num2str(psver)];
    updir = 0;
    bperpdir = dir(['..' filesep 'baselineGRID_*']);
    if isempty(bperpdir)
        bperpdir = dir(['..' filesep '..' filesep 'baselineGRID_*']);
        updir = 1;
    end

    if ~isempty(bperpdir)
        bperp_mat = zeros(n_ps, n_image - 1, 'single');
        
        loop_indices = setdiff(1:n_image, master_ix);
        
        for k = 1:length(loop_indices)
            i = loop_indices(k);
            day_str = datestr(day(i), 'yyyymmdd');
            bperp_fname = ['baselineGRID_', day_str];
            if updir, bperp_fname = ['..' filesep bperp_fname]; end
            bperp_fname = ['..' filesep bperp_fname]; 

            if ~exist([bperp_fname '.xml'], 'file')
                 fprintf('Warning: Baseline grid %s not found. Skipping.\n', bperp_fname);
                 continue;
            end
            
            bperp_grid = load_isce_internal(bperp_fname);
            [g_rows, g_cols] = size(bperp_grid);

            if g_cols == width && g_rows == len
                idx_sub = sub2ind(size(bperp_grid), ij(:,3) + 1, ij(:,2) + 1);
                bp0_ps = bperp_grid(idx_sub);
            else
                xv = linspace(0, width, g_cols);
                yv = linspace(0, len, g_rows);
                bp0_ps = interp2(xv, yv, bperp_grid, ij(:,3), ij(:,2), 'linear'); 
            end
            
            bperp_mat(:, k) = bp0_ps; 
        end
    else
        % Fallback
        bperp_vals = bperp([1:master_ix-1, master_ix+1:end]);
        bperp_mat = repmat(single(bperp_vals'), n_ps, 1);
    end

    stamps_save(bpsavename, bperp_mat);
    logit(1);

end

%% Internal Helper Functions ---
function full_path = conditional_path(filename)
    % path checking
    if exist(filename, 'file')
        full_path = filename;
    else
        full_path = ['../', filename];
    end
end

function found_path = search_file_upwards(filename)
    % search file upwards
    if exist(filename, 'file')
        found_path = filename;
        return;
    end
    up1 = ['..', filesep, filename];
    if exist(up1, 'file')
        found_path = up1;
        return;
    end
    up2 = ['..', filesep, '..', filesep, filename];
    if exist(up2, 'file')
        found_path = up2;
        return;
    end
    found_path = filename; 
end

function data_out = load_isce_internal(datafile)
% replace load_isce.m & get_parm_xml.m 
    xmlfile = [datafile, '.xml'];
    if ~exist(xmlfile, 'file')
        error([xmlfile ' does not exist']);
    end

    meta = parse_isce_xml_fast(xmlfile);
    
    width = meta.width;
    nbands = meta.number_bands;
    dtype = meta.data_type;
    scheme = meta.scheme;

    switch lower(dtype)
        case {'cfloat', 'float'}
            mat_type = 'float32';
        case 'double'
            mat_type = 'double';
        case 'byte'
            mat_type = 'uint8';
        case 'short'
            mat_type = 'short';
        otherwise
            error(['Unknown data type: ' dtype]);
    end

    fid = fopen(datafile, 'r');
    if fid == -1, error(['Cannot open ' datafile]); end

    is_complex = strcmpi(dtype, 'cfloat');
    
    if is_complex
        raw = fread(fid, [2 * width, inf], mat_type);
        fclose(fid);
        data_out = complex(raw(1:2:end, :), raw(2:2:end, :));
        data_out = data_out.'; 
    else
        if strcmpi(scheme, 'BIL') && nbands > 1
             raw = fread(fid, [width * nbands, inf], mat_type);
             data_out = raw(1:width, :).'; 
        else
             raw = fread(fid, [width, inf], mat_type);
             data_out = raw.'; 
        end
        fclose(fid);
    end
end

function meta = parse_isce_xml_fast(xmlfile)
% extract parameter form xml file
    fid = fopen(xmlfile, 'r');
    content = fread(fid, '*char')';  
    fclose(fid);

    meta.width = extract_xml_val(content, 'width', true);
    meta.length = extract_xml_val(content, 'length', true);
    meta.number_bands = extract_xml_val(content, 'number_bands', true);
    meta.data_type = extract_xml_string(content, 'data_type');
    meta.scheme = extract_xml_string(content, 'scheme');
    
    if isempty(meta.number_bands), meta.number_bands = 1; end
    if isempty(meta.scheme), meta.scheme = 'BSQ'; end
end

function val = extract_xml_val(content, tag, is_num)
    % regular expression extraction
    % fit model: name="tag">...<value>THE_VALUE</value>
    pattern = ['name="' tag '".*?<value>(.*?)</value>'];
    tokens = regexp(content, pattern, 'tokens', 'once');
    
    if isempty(tokens)
        val = [];
    else
        str_val = strtrim(tokens{1});
        if is_num
            val = str2double(str_val);
        else
            val = str_val;
        end
    end
end

function val = extract_xml_string(content, tag)
    val = extract_xml_val(content, tag, false);
    if isempty(val), val = ''; end
end