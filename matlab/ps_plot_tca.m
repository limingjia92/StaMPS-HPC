function [aps_corr, fig_name_tca, aps_flag_out] = ps_plot_tca(aps, aps_flag_in)
%PS_PLOT_TCA Select atmospheric phase screen (APS) correction data.
%   [aps_corr, fig_name_tca, aps_flag] = ps_plot_tca(aps, aps_flag)
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
%   1. Structural Refactoring: Replaced legacy if-elseif chains with a clear, 
%      performant switch-case architecture based on string keys.
%   2. String-Driven Logic: Promotes usage of readable string flags (e.g., 'a_era5') 
%      while maintaining backward compatibility via an integer-to-string mapping layer.
%   3. Extended Model Support: Added native support for ERA5 ('a_era5', 'a_e5'), 
%      Mcandis stacking ('a_mcandis'), and custom addition models ('a_add').
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert , October 2013
%   ======================================================================

    % 1. Input Standardization: Convert all inputs to string keys
    if ischar(aps_flag_in) || isstring(aps_flag_in)
        query_key = char(aps_flag_in);
    elseif isnumeric(aps_flag_in)
        % Compatibility Layer: Convert legacy integer IDs to string keys
        query_key = numeric_to_string_map(aps_flag_in);
    else
        error('ps_plot_tca:InvalidInput', 'Input flag must be a string or integer ID.');
    end
    
    % Return the resolved string flag for external reference
    aps_flag_out = query_key; 
    
    % Initialize outputs
    aps_corr = [];
    fig_name_tca = '';

    % 2. Core Logic: Switch-Case
    switch query_key
        
        % --- Linear & Powerlaw ---
        case {'a_linear', 'a_l'} % ID: 1
            % aps topo correlated linear correction
            aps_corr = aps.ph_tropo_linear;
            fig_name_tca = '(linear)';
            
        case {'a_powerlaw', 'a_p'} % ID: 2
            % aps topo correlated powerlaw correction
            aps_corr = aps.ph_tropo_powerlaw;
            fig_name_tca = '(powerlaw)';
            
        case {'a_powerlaw-k', 'a_pk'} % ID: 11
            % Spatial maps of the coefficient relating phase and tropo for power-law
            aps_corr = aps.K_tropo_powerlaw;
            fig_name_tca = '(powerlaw - spatial K map)';
            
        case {'a_linear-man', 'a_lman'} % ID: 18
            % aps topo correlated manually estimated
            aps_corr = aps.strat_corr;
            fig_name_tca = '(linear)';
            
        % --- MERIS (m) ---
        case {'a_meris', 'a_m'} % ID: 3
            % aps topo correlated meris correction
            aps_corr = aps.ph_tropo_meris;
            fig_name_tca = '(meris)';
            
        case {'a_meris-ni', 'a_mi'} % ID: 10
            % aps topo correlated MERIS (non-interpolated)
            aps_corr = aps.ph_tropo_meris_no_interp;
            fig_name_tca = '(meris)';
            
        % --- ERA-I (e) ---
        case {'a_erai', 'a_e'} % ID: 4
            % aps topo correlated ERA-I correction
            aps_corr = aps.ph_tropo_era;
            fig_name_tca = '(era)';
            
        case {'a_erai-h', 'a_eh'} % ID: 5
            % aps hydrostatic ERA-I correction
            aps_corr = aps.ph_tropo_era_hydro;
            fig_name_tca = '(era hydro)';
            
        case {'a_erai-w', 'a_ew'} % ID: 6
            % aps wet ERA-I correction 
            aps_corr = aps.ph_tropo_era_wet;
            fig_name_tca = '(era wet)';
            
        % --- WRF (w) ---
        case {'a_wrf', 'a_w'} % ID: 7
            % aps topo correlated WRF correction
            aps_corr = aps.ph_tropo_wrf;
            fig_name_tca = '(wrf)';
            
        case {'a_wrf-h', 'a_wh'} % ID: 8
            % aps hydrostatic WRF correction
            aps_corr = aps.ph_tropo_wrf_hydro;
            fig_name_tca = '(wrf hydro)';
            
        case {'a_wrf-w', 'a_ww'} % ID: 9
            % aps wet WRF correction
            aps_corr = aps.ph_tropo_wrf_wet;
            fig_name_tca = '(wrf wet)';
            
        % --- MODIS (M) ---
        case {'a_modis', 'a_M'} % ID: 12
            % aps topo correlated modis correction
            aps_corr = aps.ph_tropo_modis;
            fig_name_tca = '(modis)';
            
        case {'a_modis-ni', 'a_MI'} % ID: 13
            % aps topo correlated modis (non-interpolated)
            aps_corr = aps.ph_tropo_modis_no_interp;
            fig_name_tca = '(modis)';
            
        case {'a_recalmodis', 'a_RM'} % ID: 19
            % aps topo correlated modis recalibrated correction
            aps_corr = aps.ph_tropo_modis_recal;
            fig_name_tca = '(modis recal)';
            
        case {'a_recalmodis-ni', 'a_RMI'} % ID: 20
            % aps topo correlated modis recalibrated (non-interpolated)
            aps_corr = aps.ph_tropo_modis_no_interp_recal;
            fig_name_tca = '(modis recal)';
            
        % --- Combined Models (MERIS/MODIS + ERA) ---
        case {'a_meris+a_erai-h', 'a_m+a_eh'} % ID: 14
            % aps topo correlated MERIS plus a hydrostatic component from ERA-I
            ix_no_meris = sum(aps.ph_tropo_meris) == 0;
            aps_corr = aps.ph_tropo_meris + aps.ph_tropo_era_hydro;
            aps_corr(:, ix_no_meris) = 0; % Masking logic
            fig_name_tca = '(meris + ERA hydro)';
            
        case {'a_meris-ni+a_erai-h', 'a_mi+a_eh'} % ID: 15
            % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from ERA-I
            aps_corr = aps.ph_tropo_meris_no_interp + aps.ph_tropo_era_hydro;
            fig_name_tca = '(meris + ERA hydro)';
            
        case {'a_modis+a_erai-h', 'a_M+a_eh'} % ID: 16
            % aps topo correlated modis plus a hydrostatic component from ERA-I
            ix_no_modis = sum(aps.ph_tropo_modis) == 0;
            aps_corr = aps.ph_tropo_modis + aps.ph_tropo_era_hydro;
            aps_corr(:, ix_no_modis) = 0; % Masking logic
            fig_name_tca = '(modis + ERA hydro)';
            
        case {'a_modis-ni+a_erai-h', 'a_MI+a_eh'} % ID: 17
            % aps topo correlated modis (non-interpolated) plus a hydrostatic component from ERA-I
            aps_corr = aps.ph_tropo_modis_no_interp + aps.ph_tropo_era_hydro;
            fig_name_tca = '(modis + ERA hydro)';
            
        case {'a_recalmodis+a_erai-h', 'a_RM+a_eh'} % ID: 21
            % aps topo correlated modis recalibrated plus a hydrostatic component from ERA-I
            ix_no_modis = sum(aps.ph_tropo_modis) == 0;
            aps_corr = aps.ph_tropo_modis_recal + aps.ph_tropo_era_hydro;
            aps_corr(:, ix_no_modis) = 0; % Masking logic
            fig_name_tca = '(modis recal + ERA hydro)';
            
        case {'a_recalmodis-ni+a_erai-h', 'a_RMI+a_eh'} % ID: 22
            % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from ERA-I
            aps_corr = aps.ph_tropo_modis_no_interp_recal + aps.ph_tropo_era_hydro;
            fig_name_tca = '(modis recal + ERA hydro)';
            
        % --- Combined Models (MERIS/MODIS + WRF) ---
        case {'a_meris+a_wrf-h', 'a_m+a_wh'} % ID: 23
            % aps topo correlated MERIS plus a hydrostatic component from WRF
            ix_no_meris = sum(aps.ph_tropo_meris) == 0;
            aps_corr = aps.ph_tropo_meris + aps.ph_tropo_wrf_hydro;
            aps_corr(:, ix_no_meris) = 0; % Masking logic
            fig_name_tca = '(meris + WRF hydro)';
            
        case {'a_meris-ni+a_wrf-h', 'a_mi+a_wh'} % ID: 24
            % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from WRF
            aps_corr = aps.ph_tropo_meris_no_interp + aps.ph_tropo_wrf_hydro;
            fig_name_tca = '(meris + WRF hydro)';
            
        case {'a_modis+a_wrf-h', 'a_M+a_wh'} % ID: 25
            % aps topo correlated modis plus a hydrostatic component from WRF
            ix_no_modis = sum(aps.ph_tropo_modis) == 0;
            aps_corr = aps.ph_tropo_modis + aps.ph_tropo_wrf_hydro;
            aps_corr(:, ix_no_modis) = 0; % Masking logic
            fig_name_tca = '(modis + WRF hydro)';
            
        case {'a_modis-ni+a_wrf-h', 'a_MI+a_wh'} % ID: 26
            % aps topo correlated modis (non-interpolated) plus a hydrostatic component from WRF
            aps_corr = aps.ph_tropo_modis_no_interp + aps.ph_tropo_wrf_hydro;
            fig_name_tca = '(modis + WRF hydro)';
            
        case {'a_recalmodis+a_wrf-h', 'a_RM+a_wh'} % ID: 27
            % aps topo correlated modis recalibrated plus a hydrostatic component from WRF
            ix_no_modis = sum(aps.ph_tropo_modis) == 0;
            aps_corr = aps.ph_tropo_modis_recal + aps.ph_tropo_wrf_hydro;
            aps_corr(:, ix_no_modis) = 0; % Masking logic
            fig_name_tca = '(modis recal + WRF hydro)';
            
        case {'a_recalmodis-ni+a_wrf-h', 'a_RMI+a_wh'} % ID: 28
            % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
            aps_corr = aps.ph_tropo_modis_no_interp_recal + aps.ph_tropo_wrf_hydro;
            fig_name_tca = '(modis recal + WRF hydro)';
            
        % --- MERRA ---
        case 'a_merra' % ID: 29
            % MERRA correction
            aps_corr = aps.ph_tropo_merra;
            fig_name_tca = '(MERRA)';
            
        case 'a_merra2' % ID: 30
            % MERRA-2 correction
            aps_corr = aps.ph_tropo_merra2;
            fig_name_tca = '(MERRA-2)';
            
        case 'a_merra-h' % ID: 31
            % MERRA hydro correction
            aps_corr = aps.ph_tropo_merra_hydro;
            fig_name_tca = '(MERRA hydro)';
            
        case 'a_merra2-h' % ID: 32
            % MERRA-2 hydro correction
            aps_corr = aps.ph_tropo_merra2_hydro;
            fig_name_tca = '(MERRA-2 hydro)';
            
        case 'a_merra-w' % ID: 33
            % MERRA wet correction
            aps_corr = aps.ph_tropo_merra_wet;
            fig_name_tca = '(MERRA wet)';
            
        case 'a_merra2-w' % ID: 34
            % MERRA-2 wet correction
            aps_corr = aps.ph_tropo_merra2_wet;
            fig_name_tca = '(MERRA-2 wet)';
            
        % --- GACOS ---
        case 'a_gacos' % ID: 35
            % GACOS correction
            aps_corr = aps.ph_tropo_gacos;
            fig_name_tca = '(GACOS)';
            
        % --- NARR ---
        case 'a_narr' % ID: 36
            % NARR correction total
            aps_corr = aps.ph_tropo_narr;
            fig_name_tca = '(NARR)';
            
        case 'a_narr-h' % ID: 37
            % NARR hydro correction
            aps_corr = aps.ph_tropo_narr_hydro;
            fig_name_tca = '(NARR hydro)';
            
        case 'a_narr-w' % ID: 38
            % NARR wet correction
            aps_corr = aps.ph_tropo_narr_wet;
            fig_name_tca = '(NARR wet)';
            
        % --- ERA5 ---
        case {'a_era5', 'a_e5'} % ID: 39 
            % aps topo correlated ERA5 correction
            aps_corr = aps.ph_tropo_era5;
            fig_name_tca = '(era5)';
            
        case {'a_era5-h', 'a_e5h'} % ID: 40
            % aps hydrostatic ERA5 correction
            aps_corr = aps.ph_tropo_era5_hydro;
            fig_name_tca = '(era hydro5)';
            
        case {'a_era5-w', 'a_e5w'} % ID: 41
            % aps wet ERA5 correction (wet)
            aps_corr = aps.ph_tropo_era5_wet;
            fig_name_tca = '(era wet5)';
            
        % --- Custom / Other ---
        case {'a_mcandis', 'a_mc'} % ID: 42
            % aps mcandis stacking correction
            aps_corr = aps.aps_out;
            fig_name_tca = '(mcandis)';
            
        case 'a_add' % ID: 43
            % aps topo correlated ERA5 corection add mcandis stacking corection
            aps_corr = aps.aps_add;
            fig_name_tca = '(era5 add mcandis)';
            
        case 'none'
            % [StaMPS-HPC Safety Guard]: Missing APS flag interception
            error('StaMPS_HPC:MissingAPSFlag', ...
                  ['Atmospheric correction (''a'') requested, but no APS model specified. ', ...
                   'Please append a valid APS flag to your command (e.g., ''a_era5'', ''a_linear'').']);
            
        otherwise
            error('StaMPS_HPC:InvalidAPSFlag', 'Not a valid APS option: %s', query_key);
    end
end

% =========================================================
% Helper Function: Map Numeric ID to String Key
% =========================================================
function str_key = numeric_to_string_map(id)
    switch id
        case 1,  str_key = 'a_linear';
        case 2,  str_key = 'a_powerlaw';
        case 3,  str_key = 'a_meris';
        case 4,  str_key = 'a_erai';
        case 5,  str_key = 'a_erai-h';
        case 6,  str_key = 'a_erai-w';
        case 7,  str_key = 'a_wrf';
        case 8,  str_key = 'a_wrf-h';
        case 9,  str_key = 'a_wrf-w';
        case 10, str_key = 'a_meris-ni';
        case 11, str_key = 'a_powerlaw-k';
        case 12, str_key = 'a_modis';
        case 13, str_key = 'a_modis-ni';
        case 14, str_key = 'a_meris+a_erai-h';
        case 15, str_key = 'a_meris-ni+a_erai-h';
        case 16, str_key = 'a_modis+a_erai-h';
        case 17, str_key = 'a_modis-ni+a_erai-h';
        case 18, str_key = 'a_linear-man';
        case 19, str_key = 'a_recalmodis';
        case 20, str_key = 'a_recalmodis-ni';
        case 21, str_key = 'a_recalmodis+a_erai-h';
        case 22, str_key = 'a_recalmodis-ni+a_erai-h';
        case 23, str_key = 'a_meris+a_wrf-h';
        case 24, str_key = 'a_meris-ni+a_wrf-h';
        case 25, str_key = 'a_modis+a_wrf-h';
        case 26, str_key = 'a_modis-ni+a_wrf-h';
        case 27, str_key = 'a_recalmodis+a_wrf-h';
        case 28, str_key = 'a_recalmodis-ni+a_wrf-h';
        case 29, str_key = 'a_merra';
        case 30, str_key = 'a_merra2';
        case 31, str_key = 'a_merra-h';
        case 32, str_key = 'a_merra2-h';
        case 33, str_key = 'a_merra-w';
        case 34, str_key = 'a_merra2-w';
        case 35, str_key = 'a_gacos';
        case 36, str_key = 'a_narr';
        case 37, str_key = 'a_narr-h';
        case 38, str_key = 'a_narr-w';
        case 39, str_key = 'a_era5';
        case 40, str_key = 'a_era5-h';
        case 41, str_key = 'a_era5-w';
        case 42, str_key = 'a_mcandis';
        case 43, str_key = 'a_add';
        otherwise
            error('not a valid APS option');
    end
end