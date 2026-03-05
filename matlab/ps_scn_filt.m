function []=ps_scn_filt()
%PS_SCN_FILT Estimate spatially correlated noise in unwrapped phase
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
%   1. Topology Generation: Replaced legacy 'triangle' C-call and disk I/O 
%      with native 'delaunayTriangulation', drastically improving speed.
%   2. Vectorized Filtering & Solving: Used multi-RHS sparse matrix solving 
%      and vectorized time-domain filtering to eliminate nested loops.
%   3. KD-Tree Spatial Search: Replaced sliding window with 'rangesearch' 
%      to boost speed and fix a heuristic bug causing missing neighbors.
%   4. Parallel Processing: Implemented automatic parallel pool management 
%      and 'parfor' to leverage multi-core CPUs for spatial filtering.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

logit;
fprintf('Estimating other spatially-correlated noise (HPC Optimized)...\n')

time_win=getparm('scn_time_win',1);
deramp_ifg=getparm('scn_deramp_ifg',1);
scn_wavelength=getparm('scn_wavelength',1);
drop_ifg_index=getparm('drop_ifg_index',1);
small_baseline_flag=getparm('small_baseline_flag',1);

load psver
psname=['ps',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
scnname=['scn',num2str(psver)]; 

ps=load(psname);
uw=load(phuwname);

if strcmpi(small_baseline_flag,'y')
    unwrap_ifg_index=[1:ps.n_image];
else
    unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
end

day=ps.day(unwrap_ifg_index);
master_ix=sum(ps.master_day>ps.day)+1;
n_ifg=length(unwrap_ifg_index);
n_ps=ps.n_ps;

ph_all=single(uw.ph_uw(:,unwrap_ifg_index));
clear uw

if exist([sclaname,'.mat'],'file')
    scla=load(sclaname);
    ph_all=ph_all-single(scla.ph_scla(:,unwrap_ifg_index));
    ph_all=ph_all-repmat(single(scla.C_ps_uw),1,length(unwrap_ifg_index));
    if ~isempty(scla.ph_ramp)
        ph_all=ph_all-single(scla.ph_ramp(:,unwrap_ifg_index));
    end
    clear scla
end
ph_all(isnan(ph_all))=0;  

fprintf('   Number of points per ifg: %d\n',n_ps)

% =========================================================================
% OPTIMIZATION 1: Native MATLAB Delaunay Triangulation
% =========================================================================
fprintf('   Creating Delaunay triangulation in memory...\n')
DT = delaunayTriangulation(double(ps.xy(:,2)), double(ps.xy(:,3)));
E = edges(DT);
N = size(E, 1);
edges_nz = zeros(N, 4);
edges_nz(:, 2:3) = E;
edges_nz(:, 1) = 1:N;
disp([num2str(N),' edges created successfully.'])

%%% Deramp end interferograms
if strcmpi(deramp_ifg,'all')
    deramp_ifg=1:ps.n_ifg;
end
deramp_ifg=intersect(deramp_ifg,unwrap_ifg_index);
deramp_ix=zeros(size(deramp_ifg));
ph_ramp=zeros(n_ps,length(deramp_ifg));

if ~isempty(deramp_ifg)
    fprintf('   Deramping selected interferograms...\n')
    G=double([ones(n_ps,1),ps.xy(:,2),ps.xy(:,3)]);

    for i=1:length(deramp_ifg)
        i3=find(unwrap_ifg_index==deramp_ifg(i));
        deramp_ix(i)=i3;
        d=(ph_all(:,i3));
        m=G\double(d(:));
        ph_this_ramp=G*m;
        ph_all(:,i3)=ph_all(:,i3)-ph_this_ramp; 
        ph_ramp(:,i)=ph_this_ramp;
    end
    save(scnname,'ph_ramp')    
end

% =========================================================================
% OPTIMIZATION 2: Vectorized Time Domain Low-Pass Filtering
% =========================================================================
dph_hpt = ph_all(edges_nz(:,3),:) - ph_all(edges_nz(:,2),:);
clear ph_all
dph_lpt = zeros(size(dph_hpt), 'single');
n_edges = size(dph_hpt,1);

fprintf('   Low-pass filtering pixel-pairs in time...\n')
for i1=1:n_ifg
    time_diff_sq = (day(i1)-day)'.^2;
    weight_factor = exp(-time_diff_sq/2/time_win^2);
    weight_factor(master_ix) = 0; 
    weight_factor = weight_factor/sum(weight_factor);
    
    dph_lpt(:,i1) = single(double(dph_hpt) * weight_factor');
end

dph_hpt = dph_hpt - dph_lpt;  
clear dph_lpt

ref_ix = 1;
A = sparse([[1:n_edges]';[1:n_edges]'],[edges_nz(:,2);edges_nz(:,3)],[-ones(n_edges,1);ones(n_edges,1)]);
A = double(A(:,[1:ref_ix-1,ref_ix+1:n_ps]));
clear edges_nz

% =========================================================================
% OPTIMIZATION 3: Multi-RHS Sparse Matrix Solving
% =========================================================================
fprintf('   Solving for high-frequency (in time) pixel phase...\n')
ph_hpt_temp = A \ double(dph_hpt);
clear A dph_hpt

ph_hpt = [ph_hpt_temp(1:ref_ix-1, :); zeros(1,n_ifg); ph_hpt_temp(ref_ix:end, :)]; 
clear ph_hpt_temp

ph_hpt(:,deramp_ix) = ph_hpt(:,deramp_ix) + ph_ramp;
ph_hpt = single(ph_hpt);

% =========================================================================
% OPTIMIZATION 4: KD-Tree based Parallel Spatial Filtering
% =========================================================================
sigma_sq_times_2 = 2*scn_wavelength.^2;
patch_dist = scn_wavelength*4;
ps_xy_coords = double(ps.xy(:, 2:3)); 

fprintf('   Low-pass filtering in space (KD-Tree search)...\n')
idx_cell = rangesearch(ps_xy_coords, ps_xy_coords, patch_dist);

ph_scn = nan(n_ps, n_ifg, 'single');

% Start Parallel Pool
fprintf('   Starting parallel pool...\n')
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool;
end

parfor i = 1:n_ps
    in_range_ix = idx_cell{i};    
    xy_near = ps_xy_coords(in_range_ix, :);
    
    dist_sq = (xy_near(:,1) - ps_xy_coords(i,1)).^2 + (xy_near(:,2) - ps_xy_coords(i,2)).^2;
    
    weight_factor = exp(-dist_sq / sigma_sq_times_2);
    weight_factor = weight_factor / sum(weight_factor);
    
    ph_scn(i,:) = weight_factor' * double(ph_hpt(in_range_ix, :)); 
end

% Close Parallel Pool
fprintf('   Shutting down parallel pool...\n')
delete(poolobj);

fprintf('   Spatial filtering completed.\n')

% Re-reference to 1st PS
ph_scn = ph_scn - repmat(ph_scn(1,:), n_ps, 1); 

temp = matfile(phuwname);
ph_scn_slave = zeros(size(temp, 'ph_uw'), 'single');
ph_scn_slave(:,unwrap_ifg_index) = ph_scn;
clear ph_scn

ph_scn_slave(:,master_ix) = 0;

stamps_save(scnname, ph_scn_slave, ph_hpt, ph_ramp) 
logit(1);

end