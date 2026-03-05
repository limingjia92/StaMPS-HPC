function []=sb_invert_iono(iono_flag)
%SB_INVERT_IONO Invert ionospheric phase screen correction data of short baseline ifgs 
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
%   1. Input Standardization: Added string-based IONO flag support while 
%      preserving backward compatibility with legacy numeric IDs.
%   2. OOM Prevention: Replaced double-transpose 'lscov' with a projection 
%      operator (H_op) for massive memory reduction and speed gains.
%   3. Architecture Alignment: Refactored legacy if-condition into a 
%      clean string-driven switch-case to match APS inversion logic and 
%      support future extensibility.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: David Bekaert, November 2014 - University of Leeds
%   ======================================================================

logit;

load psver
psname=['./ps',num2str(psver)];
iononame = ['./ica' num2str(psver) '.mat'];
ionosbname = ['./ica_sb' num2str(psver) '.mat'];

if nargin<1
   iono_flag = []; 
end
ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
unwrap_ifg_index=setdiff(1:ps.n_ifg, drop_ifg_index);

% Build design matrix
G=zeros(ps.n_ifg,ps.n_image);
for i=1:ps.n_ifg
    G(i,ps.ifgday_ix(i,1))=-1;
    G(i,ps.ifgday_ix(i,2))=1;
end

if sum(abs(G(:,ps.master_ix)))==0 
   error('Apparently none of the unwrapped interferograms include the original master image')
else
    % Take out master as ref by setting to zero
    G(:,ps.master_ix)=0; 
end

% Extract non-zero columns
G2=G(unwrap_ifg_index,:);
nzc_ix=sum(abs(G2))~=0; % index for non-zero columns
G2=G2(:,nzc_ix);

if rank(G2)<size(G2,2) 
    error('There are isolated subsets (cannot be inverted w.r.t. master)')
end

if ~isempty(iono_flag)
    % Input Standardization
    if isnumeric(iono_flag)
        query_key = numeric_to_string_map_iono(iono_flag);
    else
        query_key = char(iono_flag);
    end
    
    aps=load(ionosbname);
    [aps_corr_sb, ~] = ps_plot_ica(aps, query_key);
    
    % Identify ifgs with valid APS estimates
    ix = find(nanmean(aps_corr_sb,1)==0);
    unwrap_ifg_index_new = setdiff(unwrap_ifg_index,ix);

    G3=G(unwrap_ifg_index_new,:);
    nzc_ix=sum(abs(G3))~=0; 
    G3=G3(:,nzc_ix);
    
    if rank(G3)<size(G3,2) 
        error('Isolated subsets that do not have an APS estimate detected (cannot be inverted w.r.t. master)')
    end

    aps_corr_sb= aps_corr_sb(:,unwrap_ifg_index_new);
    G2 = G3;
    
    % Verify master image has an estimated delay
    ifgday_new = ps.ifgday(unwrap_ifg_index_new,:);
    ifgs_days = unique(reshape(ifgday_new((sum(aps_corr_sb)~=0),:),[],1));

    if sum(ps.master_day==ifgs_days)~=1
        error('Master does not have an APS estimated, inversion not possible')
    end

    aps_corr=zeros(ps.n_ps,ps.n_image,'single');
    
    % Projection operator for memory reduction and speedup
    H_op = (G2' * G2) \ G2'; 
    aps_corr(:, nzc_ix) = single(double(aps_corr_sb) * H_op');
    
    % Route and save based on selected iono method
    switch query_key
        case {'i_azshift', 'i_az'} 
            ph_iono_azshift = aps_corr;       
            stamps_save(iononame, ph_iono_azshift);
        case {'i_split', 'i_s'}
            ph_iono_split = aps_corr;
            stamps_save(iononame, ph_iono_split);
        case {'i_tec', 'i_t'}
            ph_iono_tec = aps_corr;
            stamps_save(iononame, ph_iono_tec);
        otherwise
            error('sb_invert_iono: Unsupported option or mapping failure.');
    end
end

logit(1);
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
