function ps_gacos2tca(gacos_dir)
% PS_GACOS2TCA (HPC Optimized Version)
%   Calculate tropospheric delay phase for PS pixels using GACOS products.
%   Supports both PS and SBAS modes, automatically handling date formatting
%   and reference point alignment.
%
%   Usage:         ps_gacos2tca(gacos_dir)
%   Input:         gacos_dir - String, path to GACOS data (.tif or .ztd)
%   Output:        Saves 'tca2.mat' with 'ph_tropo_gacos' variable
%
%   Author:        Mingjia Li
%   Date:          April 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
% =========================================================================

    if nargin < 1 || isempty(gacos_dir)
        error('Input error: Please specify the GACOS data directory.');
    end

    %% 1. Initialization and Parameter Loading
    fprintf('Initializing GACOS processing...\n');
    
    lambda   = getparm('lambda');     

    load('psver.mat', 'psver');
    psname = ['ps', num2str(psver)];
    laname = ['la', num2str(psver)];
    out_file = ['tca',num2str(psver),'.mat'];
       
    %% 2. Load Coordinates and Interferogram Dates
    small_baseline_flag = getparm('small_baseline_flag');

    if strcmpi(small_baseline_flag, 'y')
        fprintf('Mode: SBAS. Loading ifgday matrix...\n');
        ps_data = load(psname, 'lonlat', 'n_ifg', 'ifgday', 'll0');
        ifg_pairs = ps_data.ifgday;
        n_ifg = ps_data.n_ifg;
    else
        fprintf('Mode: PS. Building interferogram pairs...\n');
        ps_data = load(psname, 'lonlat', 'master_day', 'day', 'n_ifg', 'll0');
        % Bind all dates (including master) with the master date
        ifg_pairs = [repmat(ps_data.master_day, ps_data.n_ifg, 1), ps_data.day];
        n_ifg = ps_data.n_ifg;
    end

    % Fetch reference PS indices using StaMPS-HPC standard method
    ref_ps_idx = ps_setref(ps_data); 

    ps_lon = ps_data.lonlat(:, 1);
    ps_lat = ps_data.lonlat(:, 2);
    n_ps   = length(ps_lon);

    % Load incidence angle (fallback to 39 degrees if not found)
    try
        la_data = load(laname, 'la'); 
        inc_angle = la_data.la;
    catch
        warning('Incidence angle file (%s) not found. Using default 39 deg.', laname);
        inc_angle = ones(n_ps, 1) * (39 * pi / 180); 
    end

    ph_tropo_gacos = zeros(n_ps, n_ifg);
    fprintf('Total PS points: %d | Total interferograms: %d\n\n', n_ps, n_ifg);

    %% 3. Compute Tropospheric Delay Phase
    for i = 1:n_ifg
        % Format dates to YYYYMMDD string
        if ifg_pairs(i, 1) < 1000000
            master_str = datestr(ifg_pairs(i, 1), 'yyyymmdd');
        else
            master_str = num2str(ifg_pairs(i, 1));
        end
        
        if ifg_pairs(i, 2) < 1000000
            slave_str = datestr(ifg_pairs(i, 2), 'yyyymmdd');
        else
            slave_str = num2str(ifg_pairs(i, 2));
        end

        fprintf('[%02d/%02d] Processing pair: %s - %s\n', i, n_ifg, master_str, slave_str);
        
        % Skip if master and slave are the same (zero-phase reference column)
        if strcmp(master_str, slave_str)
            continue;
        end
        
        % Read GACOS data
        [ztd_master, lon_m, lat_m] = read_gacos_data(gacos_dir, master_str);
        [ztd_slave, lon_s, lat_s]  = read_gacos_data(gacos_dir, slave_str);
        
        % Interpolate Zenith Total Delay (ZTD) to PS coordinates
        ztd_master_ps = interp2(lon_m, lat_m, ztd_master, ps_lon, ps_lat, 'linear');
        ztd_slave_ps  = interp2(lon_s, lat_s, ztd_slave, ps_lon, ps_lat, 'linear');
        
        % Core Physics Computation
        dztd_ps = ztd_slave_ps - ztd_master_ps;
        
        % Spatial differencing relative to reference point(s)
        ref_delay = mean(dztd_ps(ref_ps_idx), 'omitnan');
        dztd_ps   = dztd_ps - ref_delay;
        
        % Project ZTD to LOS and convert to phase (radians)
        dlos_ps = dztd_ps ./ cos(inc_angle);              
        ph_tropo_gacos(:, i) = dlos_ps * (4 * pi / lambda);    
    end

    %% 4. Save Output
    save(out_file, 'ph_tropo_gacos');
    fprintf('\nSuccessfully saved GACOS tropospheric phase to: %s\n', out_file);

end

% =========================================================================
% Helper Function: Read GACOS Data (.tif or .ztd/.rsc)
% =========================================================================
function [ztd_grid, lon_grid, lat_grid] = read_gacos_data(data_dir, date_str)
    tif_file = fullfile(data_dir, [date_str, '.ztd.tif']);
    bin_file = fullfile(data_dir, [date_str, '.ztd']);
    rsc_file = fullfile(data_dir, [date_str, '.ztd.rsc']);
    
    if exist(tif_file, 'file')
        [ztd_grid, R] = readgeoraster(tif_file);
        
        % Generate monotonically increasing grids for interp2
        lon_vec = R.LongitudeLimits(1) + R.CellExtentInLongitude/2 : R.CellExtentInLongitude : R.LongitudeLimits(2);
        lat_vec = R.LatitudeLimits(2) - R.CellExtentInLatitude/2 : -R.CellExtentInLatitude : R.LatitudeLimits(1);
        
        lat_vec = fliplr(lat_vec); 
        ztd_grid = flipud(ztd_grid);
        [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
        
    elseif exist(bin_file, 'file') && exist(rsc_file, 'file')
        [width, len, xfirst, yfirst, xstep, ystep] = read_rsc(rsc_file);
        fid = fopen(bin_file, 'rb', 'ieee-le');
        if fid == -1
            error('Cannot open file: %s', bin_file); 
        end
        raw_data = fread(fid, [width, len], 'float32');
        fclose(fid);
        ztd_grid = raw_data'; 
        
        lon_vec = xfirst : xstep : (xfirst + (width-1)*xstep);
        lat_vec = yfirst : ystep : (yfirst + (len-1)*ystep); 
        
        lat_vec = fliplr(lat_vec);
        ztd_grid = flipud(ztd_grid);
        [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    else
        error('GACOS data not found for date %s in %s', date_str, data_dir);
    end
    
    % Treat exactly 0 as NaN to prevent interpolation artifacts over nodata regions
    ztd_grid(ztd_grid == 0) = NaN;
end

% =========================================================================
% Helper Function: Parse .rsc header file
% =========================================================================
function [width, len, xfirst, yfirst, xstep, ystep] = read_rsc(rsc_file)
    fid = fopen(rsc_file, 'r');
    data = textscan(fid, '%s %s');
    fclose(fid);
    keys = data{1}; 
    vals = data{2};
    
    width  = str2double(vals{strcmp(keys, 'WIDTH')});
    len    = str2double(vals{strcmp(keys, 'FILE_LENGTH')});
    xfirst = str2double(vals{strcmp(keys, 'X_FIRST')});
    yfirst = str2double(vals{strcmp(keys, 'Y_FIRST')});
    xstep  = str2double(vals{strcmp(keys, 'X_STEP')});
    ystep  = str2double(vals{strcmp(keys, 'Y_STEP')});
end