function [oscilatior_corr_ifgs, oscilatior_corr_velocity] = env_oscilator_corr(ps, forced_sm_flag)
%ENV_OSCILATOR_CORR Perform oscillator drift correction for Envisat ASAR interferograms.
%
%   This function performs oscillator drift correction for Envisat interferograms 
%   based on Petar Marinkovic's presentation at ESA Living Planet 2013 in Edinburgh.
%   
%   Approximation formula:
%   dR/year = c/2 * (slantRangeTime_FAR - slantRangeTime_NEAR) * corrPerYear
%
%   'Apparent displacement' correction:
%   dR/year ~ (7.8m * 5000)*3.87e-7 ~ 0.01482m
%   Correction uses the pixel information in range and, in addition, the temporal 
%   baseline information for interferograms.
%
%   Correction is defined such that:
%   Corrected interferogram/velocity = original interferogram/velocity - correction
%
%   INPUTS:
%   ps                          Pre-loaded ps.mat structure containing ij, day, ifgday, etc.
%   forced_sm_flag              (Optional) 1 to force single-master processing, 0 otherwise.
%
%   OUTPUTS:
%   oscilatior_corr_ifgs        Correction for individual interferograms in rad.
%   oscilatior_corr_velocity    Correction in mm for the velocity.
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
%   1. Direct Struct Passing: Deprecated 'envisat_flag' and secondary I/O loads by 
%      directly passing the pre-loaded 'ps' struct, preventing memory bloat.
%   2. Simplified Platform Check: Streamlined platform validation strictly through 
%      getparm('platform'), removing legacy 'master.res' text-search logic.
%   3. Memory Prevention: Leveraged implicit expansion to return scalar 0s instead 
%      of full-sized zero matrices for non-Envisat sensors, saving substantial memory.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, 2014 - University of Leeds
%   Cite as: P. Marinkovic (PPO.labs) and Y. Larsen (NORUT)
%            Consequences of Long-Term ASAR Local Oscillator Frequency Decay - 
%            an Empirical Study of 10 Years of Data. ESA Living Planet Symposium (2013)
%   ======================================================================

if nargin < 2 
    forced_sm_flag = 0;
end

small_baseline_flag = getparm('small_baseline_flag');
if forced_sm_flag == 1
    small_baseline_flag = 'n';
end

% Check platform type
platform = getparm('platform');
if isempty(platform)
    platform = 'UNKNOWN';
end

if strcmpi(platform, 'ENVISAT') || strcmpi(platform, 'ENV')
    is_envisat = true;
    fprintf('This is Envisat, oscillator drift is being removed...\n');
else
    is_envisat = false;
end

% Compute oscillator drift correction
if is_envisat
    lambda = getparm('lambda');

    % Velocity map correction
    envisat_resolution = 7.8;                   % ground range resolution [m]
    Oscilator_drift_corr_year = 3.87e-7;        % drift correction in range [1/year]      

    % Velocity correction in mm
    oscilatior_corr_velocity = (envisat_resolution * ps.ij(:,3)) * Oscilator_drift_corr_year * 1000;

    % Interferogram correction in rad
    if strcmpi(small_baseline_flag, 'y')
        n_ifg = ps.n_ifg;
        delta_year = (ps.ifgday(:,2) - ps.ifgday(:,1)) ./ 365.25;
    else
        n_ifg = ps.n_image;
        delta_year = (ps.day - ps.master_day) ./ 365.25;
    end
    
    oscilatior_corr_ifgs = -4 .* pi ./ lambda .* repmat(oscilatior_corr_velocity, 1, n_ifg) ./ 1000 .* repmat(delta_year', ps.n_ps, 1);
    
else
    % Return scalar 0 for non-Envisat platforms to save memory
    oscilatior_corr_ifgs = 0;
    oscilatior_corr_velocity = 0;
end

end