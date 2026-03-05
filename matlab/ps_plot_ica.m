function [aps_corr, fig_name_tca, iono_flag_out] = ps_plot_ica(aps, iono_flag_in)
%PS_PLOT_ICA Select ionospheric phase screen correction data .
%   [aps_corr, fig_name_tca, iono_flag_out] = ps_plot_ica(aps, iono_flag_in)
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
%   1. Structural Refactoring: Replaced legacy if-condition with a clean, 
%      modular, and performant switch-case architecture.
%   2. Method Expansion: Integrated and standardized support for three core 
%      ionospheric correction methods: Azimuth Shift, Split-Spectrum, and TEC maps.
%   3. String-Driven Logic: Implemented robust string key parsing while 
%      maintaining backward compatibility for legacy integer IDs.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, October 2014 - University of Leeds
%   ======================================================================

    % 1. Input Standardization: Convert all inputs to string keys
    if ischar(iono_flag_in) || isstring(iono_flag_in)
        query_key = char(iono_flag_in);
    elseif isnumeric(iono_flag_in)
        % Compatibility Layer: Convert legacy integer IDs to string keys
        query_key = numeric_to_string_map_iono(iono_flag_in);
    else
        error('Input flag must be a string or integer ID.');
    end
    
    % Return the resolved string flag for external reference
    iono_flag_out = query_key; 
    
    % Initialize outputs
    aps_corr = [];
    fig_name_tca = '';

    % 2. Core Logic: Switch-Case
    switch query_key
        
        % --- Azimuth Shift Method  ---
        case {'i_azshift', 'i_az'} % ID: 1
            % aps iono azimuth shift correction
            aps_corr = aps.ph_iono_azshift;
            fig_name_tca = ' (azimuth shift method)';
            
        % --- Split-Spectrum Method ---
        case {'i_split', 'i_s'} % ID: 2
            % aps iono split-spectrum correction
            aps_corr = aps.ph_iono_split;
            fig_name_tca = ' (split-spectrum method)';
            
        % --- TEC maps Method ---
        case {'i_tec', 'i_t'} % ID: 3
            % aps iono Tec maps correction
            aps_corr = aps.ph_iono_tec;
            fig_name_tca = ' (tec maps method)';
            
        otherwise
            error('Invalid IONO option or unsupported method.');
    end
end

% =========================================================
% Helper Function: Map Numeric ID to String Key (IONO)
% =========================================================
function str_key = numeric_to_string_map_iono(id)
    switch id
        case 1
            str_key = 'i_azshift'; 
        case 2
            str_key = 'i_split'; 
        case 3
            str_key = 'i_tec';
        otherwise
            error('Invalid IONO numeric option.');
    end
end