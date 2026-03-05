function []=sb_invert_tca(aps_flag)
%SB_INVERT_TCA Invert atmospheric phase screen (APS) correction data short baseline ifgs 
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
%   1. Input Standardization: Added string-based APS flag support while 
%      preserving backward compatibility with legacy numeric IDs.
%   2. OOM Prevention: Replaced double-transpose 'lscov' with a projection 
%      operator (H_op) for massive memory reduction and speed gains.
%   3. Architecture Alignment: Refactored legacy if-elseif chains into a 
%      clean string-driven switch-case matching 'ps_plot_tca'.
%   4. Nomenclature Consistency: Renamed 'sb_invert_aps' to 'sb_invert_tca' 
%      to align with 'tca*.mat' data and its matching function 'ps_plot_tca'.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, November 2013
%   ======================================================================

logit;

load psver
psname=['./ps',num2str(psver)];
apssbname=['./tca_sb',num2str(psver) '.mat'];
apsname=['./tca',num2str(psver) '.mat'];

if nargin<1
   aps_flag = []; 
end
ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
unwrap_ifg_index=setdiff(1:ps.n_ifg, drop_ifg_index);

ref_ps=ps_setref(ps);

% Build design matrix
G=zeros(ps.n_ifg,ps.n_image);
for i=1:ps.n_ifg
    G(i,ps.ifgday_ix(i,1))=-1;
    G(i,ps.ifgday_ix(i,2))=1;
end

if sum(abs(G(:,ps.master_ix)))==0 
   error('None of the unwrapped interferograms include the original master image.')
else
    % Remove master reference by setting to zero
    G(:,ps.master_ix)=0; 
end

% Extract non-zero columns
G2=G(unwrap_ifg_index,:);
nzc_ix=sum(abs(G2))~=0; 
G2=G2(:,nzc_ix);

if rank(G2)<size(G2,2) 
    error('Isolated subsets detected (cannot be inverted w.r.t. master).')
end

if ~isempty(aps_flag)
    % --- HPC Opt 1: Input Standardization ---
    if isnumeric(aps_flag)
        query_key = numeric_to_string_map(aps_flag);
    else
        query_key = char(aps_flag);
    end
    
    aps=load(apssbname);
    [aps_corr_sb, fig_name_tca, resolved_flag] = ps_plot_tca(aps, query_key);
    
    % Identify ifgs with valid APS estimates
    ix = find(nanmean(aps_corr_sb,1)==0);
    unwrap_ifg_index_new = setdiff(unwrap_ifg_index,ix);

    G3=G(unwrap_ifg_index_new,:);
    nzc_ix=sum(abs(G3))~=0; 
    G3=G3(:,nzc_ix);
    
    if rank(G3)<size(G3,2) 
        error('Isolated subsets without APS estimates detected.')
    end

    aps_corr_sb= aps_corr_sb(:,unwrap_ifg_index_new);
    G2 = G3;
    
    % Verify master image has an estimated delay
    ifgday_new = ps.ifgday(unwrap_ifg_index_new,:);
    ifgs_days = unique(reshape(ifgday_new((sum(aps_corr_sb)~=0),:),[],1));

    if sum(ps.master_day==ifgs_days)~=1
        error('Master lacks an estimated APS. Inversion not possible.')
    end
    
    aps_corr=zeros(ps.n_ps,ps.n_image,'single');
    
    % --- HPC Opt 2: OOM Prevention via Projection Operator ---
    H_op = (G2' * G2) \ G2'; 
    aps_corr(:, nzc_ix) = single(double(aps_corr_sb) * H_op');
    
    % --- HPC Opt 3: String-Driven Architecture Alignment ---
    switch resolved_flag
        case {'a_linear', 'a_l'} 
            ph_tropo_linear = aps_corr;       
            stamps_save(apsname, ph_tropo_linear);
        case {'a_powerlaw', 'a_p'} 
            ph_tropo_powerlaw = aps_corr;       
            stamps_save(apsname, ph_tropo_powerlaw);
        case {'a_powerlaw-k', 'a_pk'} 
            K_tropo_powerlaw = aps_corr;
            stamps_save(apsname, K_tropo_powerlaw);
        case {'a_linear-man', 'a_lman'}
            strat_corr = aps_corr;
            stamps_save(apsname, strat_corr);
            
        case {'a_meris', 'a_m'} 
            ph_tropo_meris = aps_corr;       
            stamps_save(apsname, ph_tropo_meris);
        case {'a_meris-ni', 'a_mi'} 
            ph_tropo_meris_no_interp = aps_corr;
            stamps_save(apsname, ph_tropo_meris_no_interp);
            
        case {'a_erai', 'a_e'} 
            ph_tropo_era = aps_corr;
            stamps_save(apsname, ph_tropo_era);
        case {'a_erai-h', 'a_eh'} 
            ph_tropo_era_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_era_hydro);
        case {'a_erai-w', 'a_ew'} 
            ph_tropo_era_wet = aps_corr;
            stamps_save(apsname, ph_tropo_era_wet);
            
        case {'a_wrf', 'a_w'} 
            ph_tropo_wrf = aps_corr;
            stamps_save(apsname, ph_tropo_wrf);
        case {'a_wrf-h', 'a_wh'} 
            ph_tropo_wrf_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_wrf_hydro);
        case {'a_wrf-w', 'a_ww'} 
            ph_tropo_wrf_wet = aps_corr;
            stamps_save(apsname, ph_tropo_wrf_wet);
            
        case {'a_modis', 'a_M'} 
            ph_tropo_modis = aps_corr;
            stamps_save(apsname, ph_tropo_modis);
        case {'a_modis-ni', 'a_MI'} 
            ph_tropo_modis_no_interp = aps_corr;
            stamps_save(apsname, ph_tropo_modis_no_interp);
        case {'a_recalmodis', 'a_RM'} 
            ph_tropo_modis_recal = aps_corr;
            stamps_save(apsname, ph_tropo_modis_recal);
        case {'a_recalmodis-ni', 'a_RMI'} 
            ph_tropo_modis_no_interp_recal = aps_corr;
            stamps_save(apsname, ph_tropo_modis_no_interp_recal);
            
        case 'a_merra' 
            ph_tropo_merra = aps_corr;
            stamps_save(apsname, ph_tropo_merra);
        case 'a_merra2' 
            ph_tropo_merra2 = aps_corr;
            stamps_save(apsname, ph_tropo_merra2);
        case 'a_merra-h' 
            ph_tropo_merra_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_merra_hydro);
        case 'a_merra2-h' 
            ph_tropo_merra2_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_merra2_hydro);
        case 'a_merra-w' 
            ph_tropo_merra_wet = aps_corr;
            stamps_save(apsname, ph_tropo_merra_wet);
        case 'a_merra2-w' 
            ph_tropo_merra2_wet = aps_corr;
            stamps_save(apsname, ph_tropo_merra2_wet);
            
        case 'a_gacos' 
            ph_tropo_gacos = aps_corr;
            stamps_save(apsname, ph_tropo_gacos);
            
        case 'a_narr' 
            ph_tropo_narr = aps_corr;
            stamps_save(apsname, ph_tropo_narr);
        case 'a_narr-h' 
            ph_tropo_narr_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_narr_hydro);
        case 'a_narr-w' 
            ph_tropo_narr_wet = aps_corr;
            stamps_save(apsname, ph_tropo_narr_wet);
            
        case {'a_era5', 'a_e5'} 
            ph_tropo_era5 = aps_corr;
            stamps_save(apsname, ph_tropo_era5);
        case {'a_era5-h', 'a_e5h'} 
            ph_tropo_era5_hydro = aps_corr;
            stamps_save(apsname, ph_tropo_era5_hydro);
        case {'a_era5-w', 'a_e5w'} 
            ph_tropo_era5_wet = aps_corr;
            stamps_save(apsname, ph_tropo_era5_wet);
            
        case {'a_mcandis', 'a_mc'} 
            aps_out = aps_corr;
            stamps_save(apsname, aps_out);
        case 'a_add' 
            aps_add = aps_corr;
            stamps_save(apsname, aps_add);
            
        otherwise
            error('sb_invert_aps: Unsupported option or mapping failure.');
    end
end

logit(1);
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
            error('Invalid APS numeric option.');
    end
end