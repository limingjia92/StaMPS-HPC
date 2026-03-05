function [h_fig, phase_lims, h_axes_all] = ps_plot(value_type, varargin)
%PS_PLOT High-Performance Universal Rendering Engine for StaMPS
%   [h_fig, phase_lims, h_axes_all] = ps_plot(value_type, varargin)
%
%   USAGE:
%     ps_plot(value_type)
%     ps_plot(value_type, plot_flag, phase_lims, ref_ifg, ifg_list, 'OptionName')
%
%   INPUTS:
%     value_type :  String defining the variable to plot. Can be combined 
%                   with correction flags using hyphens (e.g., 'v-dao').
%                   Capitalize to force Single-Master (SM) inversion (e.g., 'V-d').
%                   [Basic Types]:
%                     'v' / 'vsb' - Mean LOS velocity (SM / SBAS)
%                     'vs'        - Velocity standard deviation
%                     'w' / 'p'   - Wrapped phase / Filtered wrapped phase
%                     'u' / 'usb' - Unwrapped phase (SM / SBAS)
%                     'hgt'       - Topographic elevation (m)
%                     'd' / 'dsb' - Spatially correlated DEM error
%                     'm'         - Master atmospheric/orbit error
%                     'o'         - Orbital deramp error
%                     's'         - Slave specific spatial error
%                     'a' / 'asb' - Stratified topo-correlated atmosphere (TCA)
%                     'i' / 'isb' - Ionospheric delay (ICA)
%                     't' / 'tsb' - Solid Earth tide
%
%   OPTIONAL POSITIONAL ARGUMENTS:
%     plot_flag  :  Background plotting mode [default = 1]
%                    -1 = Save the extracted matrices to a .mat file
%                     0 = Black background with lon/lat axes
%                     1 = White background with lon/lat axes
%     phase_lims :  1x2 vector with colormap limits (e.g., [-10 10]) [default = auto]
%     ref_ifg    :  Interferogram to reference to (0 = master, -1 = incremental)
%     ifg_list   :  List of interferogram indices to plot [default = all]
%
%   STRING FLAGS (Passed anywhere after value_type):
%     'ts'       :  Enable time-series interactive UI extraction.
%     'a_XXXX'   :  Specify APS correction model (e.g., 'a_era5', 'a_linear').
%     'i_XXXX'   :  Specify Ionosphere correction model (e.g., 'i_as').
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          March 2026
%   Version:       1.0 (HPC Overhaul)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimizations & Architecture Enhancements:
%   1. Parameter Parsing: Replaced fragile 'eval' with 'parse_plot_prm'. Fully 
%      migrated 'aps_flag' and 'iono_flag' to robust string-driven logic.
%   2. Logic Purification: Eradicated obsolete 'aps_band', 'ext_data', and RMSE 
%      validation logic to dramatically reduce codebase bloat and dependencies.
%   3. Envisat Optimization: Streamlined 'env_oscilator_corr' execution context 
%      for unified single-master and SBAS phase compensation.
%   4. OOM Prevention: Implemented aggressive 'clear' commands within the dynamic 
%      phase assembler to minimize peak memory footprints during matrix building.
%   5. Rendering Engine: Vectorized spatial mapping and simplified viewport bounds, 
%      bypassing outdated looping routines for instantaneous visualization.
%   6. TS Data Binding: Replaced volatile workspace variables with figure-bound 
%      'appdata' for encapsulated, conflict-free Time-Series UI interaction.
%   7. UI Overhaul & Profiling: Integrated an interactive dual-axis swath profiling 
%      tool and expanded the bottom control panel to manage dynamic inputs (Radius & Bins).
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

% -----------------------------------------------------------------------
% 1. Parameter parsing and environment initialization
% -----------------------------------------------------------------------
stdargin = nargin;
if stdargin < 1
    help ps_plot_hpc; error('Not enough input arguments.');
end

% Parse arguments using independent parser
prm = parse_plot_prm(varargin);

% Unpack variables explicitly
plot_flag   = prm.plot_flag;
phase_lims  = prm.phase_lims;
ref_ifg     = prm.ref_ifg;
ifg_list    = prm.ifg_list;

ts_flag     = prm.ts_flag;
aps_flag    = prm.aps_flag;   
iono_flag   = prm.iono_flag; 

textsize = 0; units = [];

load psver
psname = ['./ps', num2str(psver)];
pmname = ['./pm', num2str(psver)];
rcname = ['./rc', num2str(psver), '.mat'];
phname = ['./ph', num2str(psver), '.mat'];
phuwname = ['./phuw', num2str(psver)];
phuwsbname = ['./phuw_sb', num2str(psver)];
phuwsbresname = ['./phuw_sb_res', num2str(psver)];
scnname = ['./scn', num2str(psver)];
ifgstdname = ['./ifgstd', num2str(psver)];
apsname = ['./tca', num2str(psver)];
apssbname = ['./tca_sb', num2str(psver)];
iononame = ['./ica', num2str(psver), '.mat'];
ionosbname = ['./ica_sb', num2str(psver), '.mat'];
tidename = ['./tide', num2str(psver)];
tidesbname = ['./tide_sb', num2str(psver)];
hgtname = ['./hgt', num2str(psver)];
sclaname = ['./scla', num2str(psver)];
sclasbname = ['./scla_sb', num2str(psver)];
sclasmoothname = ['./scla_smooth', num2str(psver)];
sclasbsmoothname = ['./scla_smooth_sb', num2str(psver)];
meanvname = ['./mv', num2str(psver)];

ps = load(psname);
day = ps.day;
master_day = ps.master_day;
lonlat = ps.lonlat;
n_ps = ps.n_ps;
n_ifg = ps.n_ifg;
master_ix = sum(day < master_day) + 1;
ref_ps = 0;    
drop_ifg_index = getparm('drop_ifg_index');
small_baseline_flag = getparm('small_baseline_flag');
scla_deramp = getparm('scla_deramp');

% Retrieve plotting bounding box 
plot_lon_rg = getparm('plot_lon_rg');
plot_lat_rg = getparm('plot_lat_rg');

% Reset to empty if default [0 0] is returned
if isequal(plot_lon_rg, [0 0]) 
    plot_lon_rg = [];
end
if isequal(plot_lat_rg, [0 0]) 
    plot_lat_rg = [];
end

% Determine interferogram indices for plotting
if strcmpi(small_baseline_flag, 'y')
    unwrap_ifg_index_sb = setdiff([1:ps.n_ifg], drop_ifg_index);
    if ischar(value_type) && value_type(1)~='w' && value_type(1)~='p' && exist([phuwname,'.mat'],'file')
        warning('off','MATLAB:load:variableNotFound');
        phuw = load(phuwname,'unwrap_ifg_index_sm');
        warning('on','MATLAB:load:variableNotFound');
        if isfield(phuw,'unwrap_ifg_index_sm')
            unwrap_ifg_index = phuw.unwrap_ifg_index_sm; 
        else
            unwrap_ifg_index = [1:ps.n_image]; 
        end
    else
        unwrap_ifg_index = unwrap_ifg_index_sb;
    end
    
    if ~ischar(value_type) && isempty(ifg_list)
        ifg_list = unwrap_ifg_index_sb;
    elseif ischar(value_type) && length(value_type)>2 && ...
           (startsWith(value_type, 'usb', 'IgnoreCase', true) || ...
            startsWith(value_type, 'rsb', 'IgnoreCase', true) || ...
            startsWith(value_type, 'asb', 'IgnoreCase', true) || ...
            startsWith(value_type, 'isb', 'IgnoreCase', true)) && isempty(ifg_list)
        ifg_list = unwrap_ifg_index_sb;
    end
else
    unwrap_ifg_index = setdiff([1:ps.n_ifg], drop_ifg_index);
end

if (value_type(1)=='u' || value_type(1)=='s' || value_type(1)=='a' || value_type(1)=='i' || value_type(1)=='w') && isempty(ifg_list)
    ifg_list = unwrap_ifg_index;
end
if ~ischar(value_type) && size(value_type,2) < length(ifg_list)
    ifg_list = [1:size(value_type,2)]'; 
end 

if ischar(value_type), units = 'rad'; end

% -----------------------------------------------------------------------
% 2. Dynamic parsing and correction control
% -----------------------------------------------------------------------
forced_sm_flag = 0;

% Enforce Single-Master (SM) extraction for specific variables
if ischar(value_type)
    if strcmpi(value_type, 'vs') || startsWith(value_type, 'vs-', 'IgnoreCase', true)
        if strcmpi(small_baseline_flag,'y') 
            fprintf('Velocity std devs calculated from small baseline interferograms\n');
        else
            forced_sm_flag = 1; 
        end
    elseif strcmpi(value_type(1),'v')
        if strcmpi(small_baseline_flag,'y') && strcmp(value_type(1),'v')
            fprintf('Velocities calculated from small baseline interferograms\n');
        else
            forced_sm_flag = 1; 
        end
    elseif strcmpi(value_type(1),'d')
        if strcmpi(small_baseline_flag,'y') && strcmp(value_type(1),'d')
            fprintf('DEM error calculated from small baseline interferograms\n');
        else
            forced_sm_flag = 1; 
        end
    end
end

value_type = lower(value_type); % convert to lowecase for unified logic

if (startsWith(value_type, 'u') || startsWith(value_type, 'a')) && ~(startsWith(value_type, 'usb') || startsWith(value_type, 'asb'))
    forced_sm_flag = 1;       
end

% Envisat oscillator drift correction
platform_name = getparm('platform');
if isempty(platform_name), platform_name = 'UNKNOWN'; end

switch upper(platform_name)
    case {'ENVISAT', 'ENV'}
        is_envisat = true;
    otherwise
        is_envisat = false;
end

if is_envisat
    if forced_sm_flag == 1
        [ph_unw_eni_osci, ~] = env_oscilator_corr(ps, forced_sm_flag);
    else
        [ph_unw_eni_osci, ~] = env_oscilator_corr(ps);
    end
else
    ph_unw_eni_osci = 0;
end

% Dynamic Network Inversion (SBAS to SM Translation)
is_hgt = startsWith(value_type, 'hgt'); 

if ~is_hgt
    requires_sm_inversion = strcmpi(small_baseline_flag, 'y') && ...
                            ~contains(value_type, 'w') && ...
                            ~contains(value_type, 'sb') && ...
                            ~startsWith(value_type, 'vsb');

    if requires_sm_inversion
        % --- Troposphere Inversion (TCA) ---
        if contains(value_type, 'a') && ~exist([apsname, '.mat'], 'file')
            sb_invert_tca(aps_flag);
        end
        % --- Ionosphere Inversion (ICA) ---
        if contains(value_type, 'i') && ~exist([iononame, '.mat'], 'file')
            sb_invert_iono(iono_flag);   % Placeholder for future implementation
        end
        % --- Solid Earth Tide Inversion ---
        if contains(value_type, 't') && ~exist([tidename, '.mat'], 'file')
            sb_invert_tide();            % Placeholder for future implementation
        end
    end
end

% -----------------------------------------------------------------------
% 3. Dynamic phase assembler
% -----------------------------------------------------------------------
if ischar(value_type)
    split_idx = strfind(value_type, '-');
    if isempty(split_idx)
        base_type = value_type; corrs = '';
    else
        base_type = value_type(1:split_idx(1)-1); corrs = value_type(split_idx(1)+1:end);
    end
    is_wrapped = startsWith(base_type, 'w') || strcmp(base_type, 'p');
    is_velocity = contains(base_type, 'v');
    fig_name_dynamic = base_type;
    
    % --- Base Data Loading ---
    switch base_type
        case 'hgt'  % Topographic elevation
            hgt = load(hgtname); ph_all = hgt.hgt; clear hgt; units='m'; 
            plot_color_scheme_old = getparm('plot_color_scheme');
            if ~startsWith(plot_color_scheme_old, 'GMT_', 'IgnoreCase', true)
                setparm('plot_color_scheme', 'GMT_globe');
            end
        case 'w'  % Wrapped phase
            if exist(rcname,'file'), rc = load(rcname); ph_base = rc.ph_rc; else, rc = load(phname); ph_base = rc.ph; end
            if ref_ifg~=0, ph_reref = rc.ph_reref; end 
            clear rc; 
            if is_envisat, ph_base = ph_base .* exp(-1j*ph_unw_eni_osci); end
        case 'p'  % Patch wrapped phase
            pm = load(pmname); ph_base = pm.ph_patch ./ abs(pm.ph_patch);
            if ps.n_ifg~=size(ph_base,2), ph_base = [ph_base(:,1:ps.master_ix-1), zeros(ps.n_ps,1), ph_base(:,ps.master_ix:end)]; end
            clear pm;
        case 'u'  % Unwrapped phase (SM)
            uw = load(phuwname); ph_base = uw.ph_uw; clear uw; 
            if is_envisat, ph_base = ph_base - ph_unw_eni_osci; end
        case 'usb'  % Unwrapped phase (SBAS)
            uw = load(phuwsbname); ph_base = uw.ph_uw; clear uw;
            if is_envisat, ph_base = ph_base - ph_unw_eni_osci; end
        case 'rsb'  % Residual unwrapped phase (SBAS)
            uw = load(phuwsbname); res = load(phuwsbresname);
            ph_base = zeros(size(uw.ph_uw)); 
            ph_base(:,unwrap_ifg_index_sb) = uw.ph_uw(:,unwrap_ifg_index_sb) - res.ph_res(:,unwrap_ifg_index_sb);
            clear uw res; textsize = 0;
        case {'d', 'dsb'}  % DEM error (SM / SBAS)
            if strcmp(base_type,'dsb') || (strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0)
                scla = load(sclasbname,'K_ps_uw'); 
            else
                scla = load(sclaname,'K_ps_uw'); 
            end
            ph_all = scla.K_ps_uw; clear scla; units = 'rad/m'; 
        case 'm'  % Master atmosphere and orbit error
            scla = load(sclaname); ph_all = scla.C_ps_uw; clear scla; 
        case 'o'  % Orbital deramp error
            if strcmpi(scla_deramp, 'n')
                error('StaMPS_HPC:DerampNotComputed', 'scla_deramp flag set to n. Set to y and rerun Step 7 before using the -o plot command.');
            end
            scla = load(sclaname); ph_all = scla.ph_ramp; clear scla; 
        case 's'  % Slave specific spatial error
            scn = load(scnname); ph_all = scn.ph_scn_slave; clear scn; 
        case {'a', 'asb'}  % Tropospheric delay / APS (SM / SBAS)
            if strcmp(base_type,'asb'), aps = load(apssbname); else, aps = load(apsname); end
            [aps_corr, fig_name_tca] = ps_plot_tca(aps, aps_flag); ph_base = aps_corr; clear aps aps_corr;
            fig_name_dynamic = [base_type fig_name_tca]; 
            if strcmpi(aps_flag, 'a_powerlaw-k') || strcmpi(aps_flag, 'a_pk')
                units = 'rad/m^{\alpha}'; 
            end
        case {'i', 'isb'}  % Ionospheric delay (SM / SBAS)
            if strcmp(base_type,'isb'), iono = load(ionosbname); else, iono = load(iononame); end
            [iono_corr, fig_name_ica] = ps_plot_ica(iono, iono_flag); ph_base = iono_corr; clear iono iono_corr;
            fig_name_dynamic = [base_type fig_name_ica];
        case {'t', 'tsb'}  % Solid Earth tide (SM / SBAS)
            if strcmp(base_type,'tsb'), tide = load(tidesbname); else, tide = load(tidename); end
            ph_all = tide.ph_tide; clear tide; 
        case {'v', 'vdrop', 'vs'}   % Mean velocity or velocity std dev (SM)
            if strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0
                uw = load(phuwsbname);
            else
                uw = load(phuwname);
            end
            ph_base = uw.ph_uw; clear uw;
            if is_envisat, ph_base = ph_base - ph_unw_eni_osci; end
        case 'vsb'   % Mean velocity (SBAS)
            phuw = load(phuwsbname); ph_base = phuw.ph_uw; clear phuw;
            if is_envisat, ph_base = ph_base - ph_unw_eni_osci; end
        case 'wa'  % Wrapped atmospheric phase
            if strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0
                aps = load(apssbname);
            else
                aps = load(apsname);
            end
            [aps_corr, fig_name_tca] = ps_plot_tca(aps, aps_flag);
            ph_base = exp(-1j*aps_corr); 
            clear aps aps_corr;
            fig_name_dynamic = ['wa' fig_name_tca];
        case 'wf'  % Filtered wrapped phase
            uw = load('uw_grid'); gridix = zeros(size(uw.nzix)); gridix(uw.nzix) = [1:uw.n_ps];
            ph_base = zeros(ps.n_ps, uw.n_ifg);
            for i = 1:ps.n_ps, ph_base(i,:) = uw.ph(gridix(uw.grid_ij(i,1),uw.grid_ij(i,2)),:); end
            clear uw;
        otherwise
            if ~strcmp(base_type, 'vs'), error('Unknown basic value type.'); end
    end

    % --- Corrections loop ---
    if exist('ph_base', 'var')
        if ~isempty(corrs), fig_name_dynamic = [fig_name_dynamic '-']; end
        % check whether load "_sb" files
        use_sb_file = (strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0) || ...
                      strcmpi(base_type, 'usb') || strcmpi(base_type, 'vsb') || ...
                      strcmpi(base_type, 'rsb') || strcmpi(base_type, 'asb') || ...
                      strcmpi(base_type, 'isb') || strcmpi(base_type, 'tsb') || ...
                      strcmpi(base_type, 'dsb');
                      
        for c_idx = 1:length(corrs)
            c = corrs(c_idx);
            switch c
                case 'd' % DEM Error
                    if use_sb_file
                        if is_wrapped, scla = load(sclasbsmoothname); else, scla = load(sclasbname); end
                    else
                        if is_wrapped, scla = load(sclasmoothname); else, scla = load(sclaname); end
                    end
                    if is_wrapped, ph_base = ph_base .* exp(-1j*scla.ph_scla); else, ph_base = ph_base - scla.ph_scla; end
                    clear scla; fig_name_dynamic = [fig_name_dynamic 'd'];
                    
                case 'a' % Troposphere (TCA)
                    if use_sb_file, aps = load(apssbname); else, aps = load(apsname); end
                    [aps_corr, fig_name_tca] = ps_plot_tca(aps, aps_flag);
                    if is_wrapped
                        ph_base = ph_base .* exp(-1j*aps_corr);
                    else
                        ph_base = ph_base - aps_corr; 
                    end
                    clear aps aps_corr; fig_name_dynamic = [fig_name_dynamic 'a' fig_name_tca];
                    
                case 'o' % Orbital Deramp
                    if is_wrapped
                        if strcmpi(scla_deramp, 'n')
                            error('StaMPS_HPC:DerampNotComputed', 'Wrapped phase deramp (-o) requires scla_deramp=y in Step 7.');
                        end
                        if use_sb_file
                            if contains(corrs, 'd'), scla = load(sclasbsmoothname); else, scla = load(sclasbname); end
                        else
                            if contains(corrs, 'd'), scla = load(sclasmoothname); else, scla = load(sclaname); end
                        end
                        ph_base = ph_base .* exp(-1j*scla.ph_ramp); clear scla;
                    else
                        ph_base = ps_deramp(ps, ph_base); 
                    end
                    fig_name_dynamic = [fig_name_dynamic 'o'];
                    
                case 's' % Slave Phase
                    if strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0
                        error('StaMPS_HPC:DimensionMismatch', 'Slave phase correction cannot be applied directly to SBAS pairs. Force SM projection (e.g., V-s).');
                    end
                    scn = load(scnname);
                    if is_wrapped, ph_base = ph_base .* exp(-1j*scn.ph_scn_slave); else, ph_base = ph_base - scn.ph_scn_slave; end
                    clear scn; fig_name_dynamic = [fig_name_dynamic 's'];
                    
                case 'm' % Master Phase
                    if use_sb_file
                        if is_wrapped, scla = load(sclasbsmoothname); else, scla = load(sclasbname); end
                    else
                        if is_wrapped, scla = load(sclasmoothname); else, scla = load(sclaname); end
                    end
                    % Implicit expansion replaces repmat
                    if is_wrapped
                        ph_base = ph_base .* exp(-1j*scla.C_ps_uw); clear scla;
                    else
                        ph_base = ph_base - scla.C_ps_uw; clear scla;
                    end
                    fig_name_dynamic = [fig_name_dynamic 'm'];
                    
                case 'i' % Ionosphere (ICA)
                    if use_sb_file, iono = load(ionosbname); else, iono = load(iononame); end
                    [iono_corr, fig_name_ica] = ps_plot_ica(iono, iono_flag);
                    if is_wrapped, ph_base = ph_base .* exp(-1j*iono_corr); else, ph_base = ph_base - iono_corr; end
                    clear iono iono_corr; fig_name_dynamic = [fig_name_dynamic 'i' fig_name_ica];
                    
                case 't' % Tide
                    if use_sb_file, tide = load(tidesbname); else, tide = load(tidename); end
                    if is_wrapped, ph_base = ph_base .* exp(-1j*tide.ph_tide); else, ph_base = ph_base - tide.ph_tide; end
                    clear tide; fig_name_dynamic = [fig_name_dynamic 't'];
            end
        end
        
        % Final Formatting based on type
        if is_velocity
            ph_uw = ph_base; clear ph_base;
        else
            ph_all = ph_base; clear ph_base;
            
            % Identify absolute physical quantities to not use ps reference
            is_absolute_val = strcmpi(base_type, 'hgt') || ...
                              strcmpi(base_type, 'vs') || ...
                              (startsWith(base_type, 'a') && (strcmpi(aps_flag, 'a_pk') || strcmpi(aps_flag, 'a_powerlaw-k')));
            
            if ~is_wrapped && ~is_absolute_val
                ref_ps = ps_setref(ps);
            end
            
            if is_wrapped && exist('ph_reref','var') && ref_ifg~=0
                ph_all = ph_all .* conj(ph_reref(:,ref_ifg)); % Implicit expansion
            end
            if is_wrapped && strcmp(base_type,'wa'), ph_all(:, sum(ph_all~=1)==0) = complex(1,1); end
            if is_wrapped && contains(corrs, 'm'), ph_all(:, master_ix) = 1; end
            
            if ~is_wrapped && ~contains(value_type, 'sb'), ph_all(:, master_ix) = 0; end
            if strcmpi(base_type,'u') || strcmpi(base_type,'usb'), units = 'rad'; textsize = 0; end
        end
    end

    fig_name = fig_name_dynamic;

    % --- Velocity specific logic ---
    if is_velocity
        ref_ps = ps_setref(ps);
        is_sb_velocity = strcmpi(base_type, 'vsb') || (strcmpi(small_baseline_flag, 'y') && forced_sm_flag == 0);
        
        if ~is_sb_velocity
            % --- SM Velocity Network Inversion ---
            if ts_flag == 1 
                if exist([sclaname '.mat'],'file') == 2
                    scla = load(sclaname,'C_ps_uw');
                    if isfield(scla,'C_ps_uw')
                        ph_uw = ph_uw - scla.C_ps_uw; clear scla; 
                    end
                end
            end
            unwrap_ifg_index = setdiff(unwrap_ifg_index, ps.master_ix);
            if ~isempty(ifg_list)
                unwrap_ifg_index = intersect(unwrap_ifg_index, ifg_list); ifg_list = []; 
            end
            
            ph_uw = ph_uw(:, unwrap_ifg_index); day = day(unwrap_ifg_index,:);
            ph_uw = ph_uw - mean(ph_uw(ref_ps,:), 1, 'omitnan'); 
            
            if exist([ifgstdname,'.mat'],'file')
                ifgstd = load(ifgstdname);
                if isfield(ifgstd,'ifg_std')
                    ifgvar = (ifgstd.ifg_std*pi/181).^2; sm_cov = diag(ifgvar(unwrap_ifg_index)); 
                else
                    sm_cov = eye(length(unwrap_ifg_index)); 
                end
            else
                sm_cov = eye(length(unwrap_ifg_index));
            end
            
            G = [ones(size(day)), day - master_day]; lambda = getparm('lambda');
            if strcmp(base_type, 'vdrop')
                ph_all = zeros(size(ph_uw)); n = size(ph_uw,2);
                for i = 1:n
                    m = lscov(G([1:i-1,i+1:end],:), double(ph_uw(:,[1:i-1,i+1:n]))', sm_cov([1:i-1,i+1:end],[1:i-1,i+1:end]));
                    ph_all(:,i) = -m(2,:)' * 365.25 / 4 / pi * lambda * 1000; 
                end
            else     
                m = lscov(G, double(ph_uw'), sm_cov); ph_all = -m(2,:)' * 365.25 / 4 / pi * lambda * 1000; 
            end
            
        else
            % --- SBAS Velocity Network Inversion ---
            if ~isempty(ifg_list)
                unwrap_ifg_index_sb = intersect(unwrap_ifg_index_sb, ifg_list); ifg_list = []; 
            end
            ph_uw = ph_uw(:, unwrap_ifg_index_sb);
            phuwres = load(phuwsbresname,'sb_cov');
            if isfield(phuwres,'sb_cov')
                sb_cov = phuwres.sb_cov(unwrap_ifg_index_sb, unwrap_ifg_index_sb); 
            else
                sb_cov = eye(length(unwrap_ifg_index_sb)); 
            end
            ifgday_ix = ps.ifgday_ix(unwrap_ifg_index_sb,:);
            
            if sum(sum(isnan(ph_uw))) > 0
                for kk = 1:size(ph_uw,2)
                    ph_uw(:,kk) = ph_uw(:,kk) - mean(ph_uw(~isnan(ph_uw(:,kk)),kk), 1); 
                end
            else
                ph_uw = ph_uw - mean(ph_uw(ref_ps,:), 1);
            end
            G = [ones(size(ifgday_ix(:,1))), day(ifgday_ix(:,2)) - day(ifgday_ix(:,1))]; lambda = getparm('lambda');
            
            if contains(corrs, 'drop') || strcmp(base_type, 'vdrop') 
                ph_all = zeros(size(ph_uw)); n = size(ph_uw,2);
                for i = 1:n
                    m = lscov(G([1:i-1,i+1:end],:), double(ph_uw(:,[1:i-1,i+1:n])'), sb_cov([1:i-1,i+1:end],[1:i-1,i+1:end]));
                    ph_all(:,i) = -m(2,:)' * 365.25 / 4 / pi * lambda * 1000; 
                end
            else
                m = lscov(G, double(ph_uw'), sb_cov); ph_all = -m(2,:)' * 365.25 / 4 / pi * lambda * 1000; 
            end
        end
        
        try 
            save('mean_v.mat', 'm'); 
        catch ME
            fallback_path = '~/mean_v.mat';
            save(fallback_path, 'm');
            
            warning('StaMPS_HPC:WriteAccessDenied', ...
                    'Failed to save in current directory (%s). Velocities saved to fallback path: %s', ...
                    ME.message, fallback_path);
        end
        textsize = 0; units = 'mm/yr';

        if strcmpi(base_type, 'vs') % calculate standard deviation result using bootstrapping
            ph_uw_mm = double(ph_uw / 4 / pi * lambda * 1000)'; 
            G_yr = [G(:,1), G(:,2) / 365.25];

            if ~is_sb_velocity
                ph_all = ps_mean_v(ph_uw_mm, G_yr, sm_cov, 1500);
            else
                ph_all = ps_mean_v(ph_uw_mm, G_yr, sb_cov, 1500);
            end
            mean_v_std = ph_all;
            try 
                if exist([meanvname, '.mat'], 'file')
                    save(meanvname, 'mean_v_std', '-append');
                else
                    save(meanvname, 'mean_v_std');
                end
            catch ME
                fallback_path = '~/mean_v_std_backup.mat';
                save(fallback_path, 'mean_v_std');
                
                warning('StaMPS_HPC:WriteAccessDenied', ...
                        'Failed to save standard deviations to %s. Saved to fallback path: %s\nError: %s', ...
                        meanvname, fallback_path, ME.message);
            end
        end       

        % --- Time Series Plot Preparation ---
        if ts_flag == 1
            ts_savename = ['ps_plot_ts_', fig_name_dynamic];
            ts_matname_bound = [ts_savename, '.mat']; 
            
            if isfield(ps, 'bperp')
                bperp = ps.bperp;
            else
                ps_temp = load(psname, 'bperp');
                bperp = ps_temp.bperp;
            end
            
            ph_mm = -ph_uw * lambda * 1000 / (4*pi);
            save(ts_matname_bound, 'ph_mm', 'ph_all', 'lonlat', 'unwrap_ifg_index', 'ref_ps', ...
                 'day', 'n_ps', 'lambda', 'ifg_list', 'master_day', 'bperp');
            fprintf('Time Series matrices saved to: %s\n', ts_matname_bound);
        end
    end

else
    ph_all = value_type; fig_name = 'data'; value_type = 'data';
    if size(ph_all,2) ~= ps.n_ifg, ifg_list = []; end
end

% -----------------------------------------------------------------------
% 4. Plotting & GUI Generation
% -----------------------------------------------------------------------
if isempty(ifg_list), ifg_list = 1:size(ph_all,2); end
n_ifg_plot = length(ifg_list);

if ~isempty(plot_lon_rg), ix = lonlat(:,1)>=plot_lon_rg(1) & lonlat(:,1)<=plot_lon_rg(2); lonlat=lonlat(ix,:); end
if ~isempty(plot_lat_rg), ix = lonlat(:,2)>=plot_lat_rg(1) & lonlat(:,2)<=plot_lat_rg(2); lonlat=lonlat(ix,:); end

max_xy = llh2local([max(lonlat),0]', [min(lonlat),0]);
fig_ar = 4/3; useratio = 1; ar = max_xy(1)/max_xy(2); 

% Type setting 
n_y = ceil(sqrt((n_ifg_plot)*ar/fig_ar)); 
n_x = ceil((n_ifg_plot)/n_y); 
d_x = useratio/n_x; 
d_y = d_x/ar*fig_ar;

if d_y > useratio/n_y
    d_y = useratio/n_y; 
    d_x = d_y*ar/fig_ar; 
end
h_y = 0.95*d_y; 
h_x = h_y*ar/fig_ar;

y = 1-d_y:-d_y:0; x = 1-useratio:d_x:1-d_x; [imY,imX] = meshgrid(y,x);

if textsize == 0
    textsize = round(10*4/n_x);
    if textsize > 16, textsize = 16; elseif textsize < 8, textsize = 8; end
end
l_t = 1/9*abs(textsize)/10; h_t = 1/50*abs(textsize)/10; 

ph_disp = ph_all(:,ifg_list);
if isreal(ph_all)
    if ref_ifg ~= 0
        if ref_ifg == -1
            ph_disp = ph_disp - [ph_disp(:,1), ph_disp(:,1:end-1)]; 
        else
            ph_disp = ph_disp - ph_all(:,ref_ifg); % Implicit expansion
        end
    else
        ref_ifg = master_ix;
    end

    is_absolute_val = strcmpi(base_type, 'hgt') || strcmpi(base_type, 'vs') || ...
                     (startsWith(base_type, 'a') && (strcmpi(aps_flag, 'a_pk') || strcmpi(aps_flag, 'a_powerlaw-k')));
                     
    if any(ref_ps ~= 0) && ~is_absolute_val
        ref_ph = (ph_disp(ref_ps,:)); mean_ph = zeros(1,size(ph_disp,2));
        for i = 1:size(ph_disp,2)
            mean_ph(i) = mean(ref_ph(~isnan(ref_ph(:,i)),i), 1);
            if isnan(mean_ph(i))
                mean_ph(i) = 0; 
                fprintf('Interferogram (%d) does not have a reference area\n', i); 
            end
        end
        clear i; ph_disp = ph_disp - mean_ph;
    end

    phsort = sort(ph_disp(~isnan(ph_disp)));
    
    if isempty(phase_lims)
        if ~isempty(phsort)
            min_idx = max(1, ceil(length(phsort)*0.001));
            max_idx = min(length(phsort), round(length(phsort)*0.999));
            minph = phsort(min_idx); 
            maxph = phsort(max_idx); 
            phase_lims = [minph, maxph];
        else
            fprintf('Interferograms do not contain valid data.\n');
            phase_lims = [-pi, pi];
        end
    end
else
    if ref_ifg == 0, ref_ifg = master_ix; elseif ref_ifg == -1, ph_disp = ph_disp .* conj([ph_disp(:,1), ph_disp(:,1:end-1)]); end
    if ref_ps ~= 0
        ph_disp = ph_disp ./ abs(ph_disp); ref_ph = (ph_disp(ref_ps,:)); mean_ph = zeros(1,size(ph_disp,2));
        for i = 1:size(ph_disp,2), mean_ph(i) = sum(ref_ph(~isnan(ref_ph(:,i)),i)); end
        clear i; ph_disp = ph_disp .* conj(mean_ph);
    end
    phase_lims = [-pi, pi];
end

% --- Export and Foreground GUI Rendering ---
if plot_flag == -1 % for data file output
    h_fig = []; phase_lims = []; h_axes_all = []; savename = ['ps_plot_', value_type];
    try 
        stamps_save(savename, ph_disp, ifg_list); 
    catch
        stamps_save(['~/',savename], ph_disp, ifg_list); 
        fprintf('Warning: Read access only, values in home directory instead\n'); 
    end
else
    h_fig = figure; set(gcf,'renderer','zbuffer','name',fig_name);
    try h_fig.WindowState = 'maximized'; catch, end

    i_im = 0; 
    if size(ifg_list,1) > 1, ifg_list = ifg_list'; end
    for i = ifg_list
        i_im = i_im + 1;
        if n_ifg_plot > 1
            h_axes = axes('position',[imX(i_im),imY(i_im),h_x,h_y]);
            set(h_axes,'Xcolor',[1 1 1]); set(h_axes,'Ycolor',[1 1 1]);
            h_axes_all(i) = h_axes; clear h_axes
        else
            h_axes_all(i) = gca;
        end
        
        % main plotting fuction
        ps_plot_ifg(ps, ph_disp(:,i_im), plot_flag, phase_lims, plot_lon_rg, plot_lat_rg);

        box on
        if n_ifg_plot > 1, set(gca,'yticklabel',[]); set(gca,'xticklabel',[]); end
        
        xlim_v = get(gca,'xlim'); x_t = (h_x-l_t)/2/h_x*(xlim_v(2)-xlim_v(1))+xlim_v(1);
        ylim_v = get(gca,'ylim');
        if textsize > 0, y_t = (h_y-1.2*h_t)/h_y*(ylim_v(2)-ylim_v(1))+ylim_v(1); else, y_t = (0.5*h_t)/h_y*(ylim_v(2)-ylim_v(1))+ylim_v(1); end
        
        % Render time/ID labels
        if textsize ~= 0 && size(day,1) == size(ph_all,2) && (strcmpi(small_baseline_flag,'n') || strcmpi(value_type(1),'s'))
            t = text(x_t, y_t, datestr(day(i),'yyyymmdd')); set(t,'fontweight','bold','color',[0 0 0.004],'fontsize',abs(textsize));
        elseif textsize ~= 0 && size(day,1) == size(ph_all,2) && forced_sm_flag == 1
            t = text(x_t, y_t, datestr(day(i),'yyyymmdd')); set(t,'fontweight','bold','color',[0 0 0.004],'fontsize',abs(textsize));
        elseif textsize ~= 0 && ps.n_ifg == size(ph_all,2) && strcmpi(small_baseline_flag,'y')
            t = text(x_t, y_t, ['      ifg ' num2str(i)]); set(t,'fontweight','bold','color',[0 0 0.004],'fontsize',abs(textsize));
        end
        
        % Render Colorbar
        if i == ref_ifg || (isempty(intersect(ref_ifg,ifg_list)) && i == ifg_list(1))
            if n_ifg_plot > 1
                h = colorbar('South'); xlim_c = get(h,'xlim'); set(h,'xlim',[xlim_c(2)-64,xlim_c(2)]);
            else
                h = colorbar('peer',gca); ylim_c = get(h,'ylim'); set(h,'ylim',[ylim_c(2)-64,ylim_c(2)]);
            end
            
            if diff(phase_lims) > 1 || diff(phase_lims) == 0, plotlims = round(phase_lims*10)/10; else, limorder = ceil(-log10(diff(phase_lims)))+2; plotlims = round(phase_lims*10^limorder)/10^limorder; end
            
            if n_ifg_plot > 1
                set(h,'xtick',[xlim_c(2)-64,xlim_c(2)],'Xticklabel',plotlims,'xcolor','k','ycolor',[0 0 0.004],'fontweight','bold','color',[0 0 0.004],'FontSize',abs(textsize))
                h_xl = xlabel(h,units); pos = get(h_xl,'position'); pos(2) = pos(2)/2.2; set(h_xl,'position',pos,'FontSize',abs(textsize));
            else
                set(h,'ytick',[ylim_c(2)-64,ylim_c(2)],'yticklabel',plotlims,'xcolor','k','ycolor',[0 0 0.004],'fontweight','bold','color',[0 0 0.004],'FontSize',abs(textsize))
                set(get(h,'ylabel'),'String',units,'FontSize',abs(textsize))  
            end
        end
    end
end

if exist('plot_color_scheme_old','var') == 1, setparm('plot_color_scheme',plot_color_scheme_old); end
fprintf('Color Range: %g to %g %s\n', phase_lims(1), phase_lims(2), units);

% =======================================================================
% Time Series & Profile Interactive UI Data Binding
% =======================================================================
if ts_flag == 1
    figure(h_fig);
    
    clear ph_base ph_all ph_uw G sm_cov sb_cov ph_mm; 
    
    if exist('ts_matname_bound', 'var')
        setappdata(h_fig, 'ts_matname', ts_matname_bound);
    end
    
    fig_bg_color = get(h_fig, 'Color');
    
    % --- Group 1: Action Buttons (Left aligned) ---
    uicontrol('Style', 'pushbutton', 'String', 'TS plot', ...
              'Units', 'normalized', 'Position', [0.05, 0.02, 0.12, 0.05], ...
              'Callback', @(~,~) ts_plot(h_fig));   % TS PLOT button
              
    uicontrol('Style', 'pushbutton', 'String', 'TS double diff.', ...
              'Units', 'normalized', 'Position', [0.18, 0.02, 0.15, 0.05], ...
              'Callback', @(~,~) ts_plotdiff(h_fig));   % TS Diff PLOT button
              
    uicontrol('Style', 'pushbutton', 'String', 'Profile plot', ...
              'Units', 'normalized', 'Position', [0.34, 0.02, 0.12, 0.05], ...
              'Callback', @(~,~) profile_plot(h_fig));  % Profile PLOT button

    % --- Group 2: Parameter Settings (Right aligned) ---
    % Radius Settings, used for ts_plot/ts_plotdiff/profile_plot
    uicontrol('Style', 'text', 'String', 'Radius:', ...
              'Units', 'normalized', 'Position', [0.49, 0.015, 0.06, 0.04], ...
              'BackgroundColor', fig_bg_color, 'HorizontalAlignment', 'right', ...
              'FontSize', 10);
              
    uicontrol('Style', 'edit', 'String', '500', ...
              'Units', 'normalized', 'Position', [0.56, 0.02, 0.06, 0.05], ...
              'BackgroundColor', [1 1 1], 'Tag', 'RadiusEditBox');
              
    uicontrol('Style', 'text', 'String', 'm', ...
              'Units', 'normalized', 'Position', [0.62, 0.015, 0.02, 0.04], ...
              'BackgroundColor', fig_bg_color, 'HorizontalAlignment', 'left', ...
              'FontSize', 10);
              
    % Bins Settings, used for profile_plot
    uicontrol('Style', 'text', 'String', 'Bins:', ...
              'Units', 'normalized', 'Position', [0.66, 0.015, 0.04, 0.04], ...
              'BackgroundColor', fig_bg_color, 'HorizontalAlignment', 'right', ...
              'FontSize', 10);
              
    uicontrol('Style', 'edit', 'String', '20', ...
              'Units', 'normalized', 'Position', [0.71, 0.02, 0.05, 0.05], ...
              'BackgroundColor', [1 1 1], 'Tag', 'BinsEditBox');
end

end

% =========================================================================
% SUBFUNCTION: Input Argument Parser
% =========================================================================
function prm = parse_plot_prm(args)
%PARSE_PLOT_PRM Parse input arguments for ps_plot and return a struct

    % Initialize default parameter structure
    prm.ts_flag   = 0;
    prm.aps_flag  = 'none';  
    prm.iono_flag = 'none';  

    % Separate custom string flags from standard positional arguments
    clean_args = {};
    for k = 1:length(args)
        arg = args{k};
        
        if ischar(arg) || isstring(arg)
            arg_str = char(arg);
            
            if strcmp(arg_str, 'ts')
                prm.ts_flag = 1;
            elseif startsWith(arg_str, 'a_')
                prm.aps_flag = arg_str;
            elseif startsWith(arg_str, 'i_')
                prm.iono_flag = arg_str; 
            else
                % Keep standard parameter strings
                clean_args{end+1} = arg_str;
            end
        else
            clean_args{end+1} = arg;
        end
    end

    % Standard input parser
    p = inputParser;
    p.KeepUnmatched = true; 
    
    addOptional(p, 'plot_flag', 1, @isnumeric);
    addOptional(p, 'phase_lims', []);
    addOptional(p, 'ref_ifg', 0, @isnumeric);
    addOptional(p, 'ifg_list', []);

    parse(p, clean_args{:});

    % Map parsed results to the output struct
    fields = fieldnames(p.Results);
    for i = 1:numel(fields)
        prm.(fields{i}) = p.Results.(fields{i});
    end

    % Use isequal to prevent logical array errors
    if isequal(prm.phase_lims, 0)
        prm.phase_lims = []; 
    end
end