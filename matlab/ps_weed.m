function []=ps_weed(all_da_flag,no_weed_adjacent,no_weed_noisy)
% PS_WEED (HPC Optimized Version)
%   Weeds out neighboring PS and save those kept to new version.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          March 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   PERFORMANCE BENCHMARK (On Test Dataset):
%   - Original:    ~1218 seconds
%   - Optimized:   ~370 seconds 
%   - Speedup:     ~3.3x (Approx. 70% reduction in execution time)
%
%   OPTIMIZATION HIGHLIGHTS:
%   1. Graph Theory Adjacency: Replaced nested cell arrays and 'while' loops 
%      with MATLAB's 'graph' and 'conncomp' for instant connected component 
%      clustering.
%   2. C-MEX Offloading: Extracted the most computationally expensive loop 
%      (phase multiplication, wrapping, and 'lscov' fitting) into a single 
%      C-MEX call ('smooth_arcs_mex').
%   3. Parallelization: Integrated OpenMP multi-threading within the C-MEX 
%      function for safe, shared-memory phase filtering across arcs.
%   4. Feature Extension: Added support for weeding and saving 'head' 
%      (heading) variables for 3D deformation decomposition.
%
%   COMPILATION REQUIREMENTS:
%   mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" smooth_arcs_mex.c
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================

% ================= Check MEX =================
if exist('smooth_arcs_mex', 'file') ~= 3
    error('Missing MEX: smooth_arcs_mex.');
end
% =============================================

logit;
logit('Weeding selected pixels...')

if nargin<1
    all_da_flag=0;
end

time_win=getparm('weed_time_win',1);
weed_standard_dev=getparm('weed_standard_dev',1);
weed_max_noise=getparm('weed_max_noise',1);
weed_zero_elevation=getparm('weed_zero_elevation',1);
weed_neighbours=getparm('weed_neighbours',1);
drop_ifg_index=getparm('drop_ifg_index');
small_baseline_flag=getparm('small_baseline_flag',1);

if nargin<2
    if strcmpi(weed_neighbours,'y')
        no_weed_adjacent=0;
    else
        no_weed_adjacent=1;
    end
end
if nargin<3
    if weed_standard_dev>=pi && weed_max_noise>=pi
        no_weed_noisy=1;
    else
        no_weed_noisy=0;
    end
end

% Load previous version information
load psver
psname=['ps',num2str(psver)];
pmname=['pm',num2str(psver)];
phname=['ph',num2str(psver)];
selectname=['select',num2str(psver)];
hgtname=['hgt',num2str(psver),'.mat'];
laname=['la',num2str(psver),'.mat'];
incname=['inc',num2str(psver),'.mat'];
bpname=['bp',num2str(psver),'.mat'];
headname=['head',num2str(psver),'.mat'];

% Load other potential PS variables if all_da_flag is active
psothername='ps_other';
pmothername='pm_other';
selectothername='select_other';
hgtothername='hgt_other';
laothername='la_other';
incothername='inc_other';
bpothername='bp_other';
headothername='head_other';

ps=load(psname);
ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
sl=load(selectname);

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

day=ps.day;
bperp=ps.bperp;
master_day=ps.master_day;

% Retrieve keeping index depending on the selection step output
if isfield(sl,'keep_ix')
    ix2=sl.ix(sl.keep_ix);
    K_ps2=sl.K_ps2(sl.keep_ix);
    C_ps2=sl.C_ps2(sl.keep_ix);
    coh_ps2=sl.coh_ps2(sl.keep_ix);
else
    ix2=sl.ix2;
    K_ps2=sl.K_ps2;
    C_ps2=sl.C_ps2;
    coh_ps2=sl.coh_ps2;
end

ij2=ps.ij(ix2,:);
xy2=ps.xy(ix2,:);
ph2=ph(ix2,:);
lonlat2=ps.lonlat(ix2,:);

pm=load(pmname);
ph_patch2=pm.ph_patch(ix2,:); 
if isfield(sl,'ph_res2')
    ph_res2=sl.ph_res2(sl.keep_ix,:);
else
    ph_res2=[];
end
clear pm
clear sl

if isfield(ps,'ph')
    ps=rmfield(ps,'ph');
end
ps=rmfield(ps,{'xy','ij','lonlat','sort_ix'});

% Append high D_A points if requested
if all_da_flag~=0
    pso=load(psothername);
    slo=load(selectothername);
    ix_other=slo.ix_other;
    n_ps_other=sum(ix_other);
    K_ps_other2=pso.K_ps_other(ix_other);
    C_ps_other2=pso.C_ps_other(ix_other);
    coh_ps_other2=pso.coh_ps_other(ix_other);
    ph_res_other2=pso.ph_res_other(ix_other,:);
    ij2=[ij2;pso.ij_other(ix_other,:)];
    xy2=[xy2;pso.xy_other(ix_other,:)];
    ph2=[ph2;pso.ph_other(ix_other,:)];
    lonlat2=[lonlat2;pso.lonlat_other(ix_other,:)];
    clear pso slo

    pmo=load(pmothername);
    ph_patch_other2=pmo.ph_patch_other(ix_other,:);
    clear pm

    K_ps2=[K_ps2;K_ps_other2];
    C_ps2=[C_ps2;C_ps_other2];
    coh_ps2=[coh_ps2;coh_ps_other2];
    ph_patch2=[ph_patch2;ph_patch_other2];
    ph_res2=[ph_res2;ph_res_other2];
else
    n_ps_other=0;
end

if exist(hgtname,'file')
    ht=load(hgtname);
    hgt=ht.hgt(ix2);
    clear ht
    if all_da_flag~=0
        hto=load(hgtothername);
        hgt=[hgt;hto.hgt_other(ix_other)];
        clear hto
    end
end

n_ps_low_D_A=length(ix2);
n_ps=n_ps_low_D_A + n_ps_other;
ix_weed=logical(ones(n_ps,1));
logit([num2str(n_ps_low_D_A),' low D_A PS, ',num2str(n_ps_other),' high D_A PS']);

% =========================================================================
% FIND AND DROP ADJACENT PIXELS (Graph Theory Optimized)
% =========================================================================
if no_weed_adjacent==0
    step_name='INITIALISE NEIGHBOUR MATRIX AND SELECT BEST (GRAPH OPTIMIZED)';
    fprintf([step_name '\n'])
    
    % Shift indices to avoid negative/zero boundaries
    ij_shift=ij2(:,2:3)+repmat([2,2]-min(ij2(:,2:3)),n_ps,1);
    neigh_ix=zeros(max(ij_shift(:,1))+1,max(ij_shift(:,2))+1);
    miss_middle=logical(ones(3));
    miss_middle(2,2)=0;
    
    % Populate the neighbor matrix with point indices
    for i=1:n_ps
        neigh_this=neigh_ix(ij_shift(i,1)-1:ij_shift(i,1)+1,ij_shift(i,2)-1:ij_shift(i,2)+1);
        neigh_this(neigh_this==0&miss_middle)=i;
        neigh_ix(ij_shift(i,1)-1:ij_shift(i,1)+1,ij_shift(i,2)-1:ij_shift(i,2)+1)=neigh_this;
    end
    
    % Replace slow cell arrays and while loop with Graph Connected Components
    linear_idx = sub2ind(size(neigh_ix), ij_shift(:,1), ij_shift(:,2));
    cluster_id = neigh_ix(linear_idx);
    
    valid_edges = (cluster_id ~= 0);
    u_nodes = find(valid_edges);
    v_nodes = cluster_id(valid_edges);
    
    % Group all adjacent points instantly using graph theory
    G_adj = graph(u_nodes, v_nodes, [], n_ps);
    bins = conncomp(G_adj); 
    
    % Sort by coherence descending
    [~, sorted_idx] = sort(coh_ps2, 'descend');
    sorted_bins = bins(sorted_idx);
    
    % Keep only the point with the highest coherence in each connected component
    [~, keep_idx_in_sorted] = unique(sorted_bins, 'stable');
    keep_ps = sorted_idx(keep_idx_in_sorted);
    
    ix_weed = false(n_ps, 1);
    ix_weed(keep_ps) = true;
    
    logit([num2str(sum(ix_weed)),' PS kept after dropping adjacent pixels']);
end

% Optional weeding of sea level points
if strcmpi(weed_zero_elevation,'y') & exist('hgt','var')
    sea_ix=hgt<1e-6;
    ix_weed(sea_ix)=false;
    logit([num2str(sum(ix_weed)),' PS kept after weeding zero elevation']);
end
xy_weed=xy2(ix_weed,:);
n_ps=sum(ix_weed);

% =========================================================================
% REMOVE DUPLICATED POINTS
% =========================================================================
ix_weed_num=find(ix_weed); 
[dummy,I]=unique(xy_weed(:,2:3),'rows');
dups=setxor(I,[1:sum(ix_weed)]'); 

for i=1:length(dups)
    dups_ix_weed=find(xy_weed(:,2)==xy_weed(dups(i),2)&xy_weed(:,3)==xy_weed(dups(i),3));
    dups_ix=ix_weed_num(dups_ix_weed);
    [dummy,I]=max(coh_ps2(dups_ix));
    ix_weed(dups_ix([1:end]~=I))=0; 
end

if ~isempty(dups)
    xy_weed=xy2(ix_weed,:);
    logit(sprintf('%d PS with duplicate lon/lat dropped\n',length(dups)'))
end

n_ps=sum(ix_weed);
ix_weed2=true(n_ps,1);

% =========================================================================
% WEED NOISY PIXELS (C-MEX Accelerated)
% =========================================================================
ps_std=zeros(n_ps,1);
ps_max=zeros(n_ps,1);

if n_ps~=0
    if no_weed_noisy==0

        step_name='DROP NOISY';
        DT=delaunayTriangulation(double(xy_weed(:,2:3)));
        edgs=edges(DT);
        n_edge=size(edgs,1);

        % Pre-process phase and correct for master noise/DEM error
        ph_weed=ph2(ix_weed,:).*exp(-j*(K_ps2(ix_weed)*bperp'));  
        ph_weed=ph_weed./abs(ph_weed);
        if ~strcmpi(small_baseline_flag,'y')
            ph_weed(:,ps.master_ix)=exp(j*(C_ps2(ix_weed)));  
        end
        
        edge_std=zeros(n_edge,1);
        edge_max=zeros(n_edge,1);
        
        % Calculate phase difference across arcs
        dph_space=(ph_weed(edgs(:,2),:).*conj(ph_weed(edgs(:,1),:)));
        dph_space=dph_space(:,ifg_index);
        n_use=length(ifg_index);
        
        for i=1:length(drop_ifg_index)
            if strcmpi(small_baseline_flag,'y')
                logit(sprintf('%s-%s dropped from noise estimation',datestr(ps.ifgday(drop_ifg_index(i),2)),datestr(ps.ifgday(drop_ifg_index(i),2))));
            else
                logit(sprintf('%s dropped from noise estimation',datestr(day(drop_ifg_index(i)))));
            end
        end

        if ~strcmpi(small_baseline_flag,'y')
          logit(sprintf('Estimating noise for all arcs (C-MEX Accelerated)...'))
          
          % Global Matrix Pre-computation for spatial/temporal weighting
          W = zeros(n_use, n_use);
          for i1 = 1:n_use
              time_diff_w = (day(ifg_index(i1)) - day(ifg_index))';
              weight_w = exp(-(time_diff_w.^2)/2/time_win^2);
              W(:, i1) = weight_w / sum(weight_w);
          end
          
          W2 = W;
          W2(1:n_use+1:end) = 0; % Zero the diagonal to leave out the point itself
          
          % Calculate the global mean phase (Double precision required for C-MEX interface)
          dph_mean_all = double(dph_space * W);
          dph_smooth2 = single(dph_space * W2);
          
          % Cast variables explicitly to ensure type safety in C-MEX execution
          day_ifg_double = double(day(ifg_index));
          dph_space_double = double(dph_space);
          W_double = double(W);
          
          % Execute highly optimized C-MEX function (Replaces lscov loop)
          dph_smooth = smooth_arcs_mex(dph_space_double, dph_mean_all, day_ifg_double, W_double);
           
           dph_noise=angle(dph_space.*conj(dph_smooth));
           dph_noise2=angle(dph_space.*conj(dph_smooth2));
           ifg_var=var(dph_noise2,0,1);
           
           % Estimate and subtract arc DEM error
           K=lscov(bperp(ifg_index),double(dph_noise)',1./double(ifg_var))'; 
           dph_noise=dph_noise-K*bperp(ifg_index)';
           
           clear dph_space dph_smooth dph_smooth2 dph_noise2
           edge_std=std(dph_noise,0,2);
           edge_max=max(abs(dph_noise),[],2);
           clear dph_noise
        else
           ifg_var=var(dph_space,0,1);
           K=lscov(bperp(ifg_index),double(dph_space)',1./double(ifg_var))';
           dph_space=dph_space-K*bperp(ifg_index)';
           edge_std=std(angle(dph_space),0,2);
           edge_max=max(abs(angle(dph_space)),[],2);
           clear dph_space
        end

        logit(sprintf('Estimating max noise for all pixels...'))
        
        % Vectorized mapping of arc-level noise back to individual PS nodes
        all_nodes = [edgs(:,1); edgs(:,2)];
        all_std_vals = [edge_std; edge_std];
        all_max_vals = [edge_max; edge_max];
        
        if isa(all_std_vals, 'single')
            fill_val = single(inf);
        else
            fill_val = inf;
        end
        
        % Accumulate the minimum noise value for each node
        ps_std = accumarray(all_nodes, all_std_vals, [n_ps 1], @min, fill_val);
        ps_max = accumarray(all_nodes, all_max_vals, [n_ps 1], @min, fill_val);
        
        ps_std = single(ps_std);
        ps_max = single(ps_max);
        
        % Filter based on thresholds
        ix_weed2=ps_std<weed_standard_dev&ps_max<weed_max_noise;
        ix_weed(ix_weed)=ix_weed2;
        n_ps=sum(ix_weed);

        logit([num2str(n_ps),' PS kept after dropping noisy pixels']);
    end
end

% Keep information about the number of PS left for StaMPS logging
if exist('no_ps_info.mat','file')~=2
   stamps_step_no_ps = zeros([5 1 ]);      
else
   load('no_ps_info.mat');
   stamps_step_no_ps(4:end)=0;
end
if n_ps==0
   fprintf('***No PS points left. Updating the stamps log for this****\n')
   stamps_step_no_ps(4)=1;
end
save('no_ps_info.mat','stamps_step_no_ps')

% =========================================================================
% SAVE THE RESULTS
% =========================================================================
weedname=['weed',num2str(psver)];
stamps_save(weedname,ix_weed,ix_weed2,ps_std,ps_max,ifg_index)

coh_ps=coh_ps2(ix_weed);
K_ps=K_ps2(ix_weed);
C_ps=C_ps2(ix_weed);
ph_patch=ph_patch2(ix_weed,:);
if ~isempty(ph_res2)
    ph_res=ph_res2(ix_weed,:);
else
    ph_res=ph_res2;
end

pmname=['pm',num2str(psver+1)];
stamps_save(pmname,ph_patch,ph_res,coh_ps,K_ps,C_ps)

clear ph_patch ph_res coh_ps K_ps C_ps ph_patch2 ph_res2 coh_ps2 K_ps2 C_ps2

ph2=ph2(ix_weed,:);
ph=ph2;
phname=['ph',num2str(psver+1)];
stamps_save(phname,ph)

clear ph

xy2=xy2(ix_weed,:);
ij2=ij2(ix_weed,:);
lonlat2=lonlat2(ix_weed,:);
ps.xy=xy2;
ps.ij=ij2;
ps.lonlat=lonlat2;
ps.n_ps=size(ph2,1);
psname=['ps',num2str(psver+1)];
save(psname,'-struct','ps');
clear ps xy2 ij2 lonlat2

if exist(hgtname,'file')
    hgt=hgt(ix_weed);
    stamps_save(['hgt',num2str(psver+1),'.mat'],hgt);
    clear hgt
end

if exist(laname,'file')
    la=load(laname);
    la=[la.la(ix2)];
    if all_da_flag~=0
        lao=load(laothername);
        la=[la;lao.la_other(ix_other)];
        clear lao
    end
    la=la(ix_weed);
    stamps_save(['la',num2str(psver+1),'.mat'],la);
    clear la
end

if exist(incname,'file')
    inc=load(incname);
    inc=[inc.inc(ix2)];
    if all_da_flag~=0
        inco=load(incothername);
        inc=[inc;inco.inc_other(ix_other)];
        clear inco
    end
    inc=inc(ix_weed);
    stamps_save(['inc',num2str(psver+1),'.mat'],inc);
    clear inc
end

if exist(bpname,'file')
    bp=load(bpname);
    bperp_mat=[bp.bperp_mat(ix2,:)];
    clear bp
    if all_da_flag~=0
        bpo=load(bpothername);
        bperp_mat=[bperp_mat;bpo.bperp_other(ix_other,:)];
        clear bpo
    end
    bperp_mat=bperp_mat(ix_weed,:);
    stamps_save(['bp',num2str(psver+1),'.mat'],bperp_mat);
end

if exist(headname,'file')
    head=load(headname);
    head=[head.head(ix2)];
    if all_da_flag~=0
        heado=load(headothername);
        head=[head;heado.heading_other(ix_other)];
        clear heado
    end
    head=head(ix_weed);
    stamps_save(['head',num2str(psver+1),'.mat'],head);
    clear head
end

% Cleanup old intermediate files
if exist(['scla_smooth',num2str(psver+1),'.mat'],'file')
    delete(['scla_smooth',num2str(psver+1),'.mat'])
end
if exist(['scla',num2str(psver+1),'.mat'],'file')
    delete(['scla',num2str(psver+1),'.mat'])
end
if exist(['scla_smooth_sb',num2str(psver+1),'.mat'],'file')
    delete(['scla_smooth_sb',num2str(psver+1),'.mat'])
end
if exist(['scla_sb',num2str(psver+1),'.mat'],'file')
    delete(['scla_sb',num2str(psver+1),'.mat'])
end
if exist(['aps',num2str(psver+1),'.mat'],'file')
    delete(['aps',num2str(psver+1),'.mat'])
end
if exist(['scn',num2str(psver+1),'.mat'],'file')
    delete(['scn',num2str(psver+1),'.mat'])
end

setpsver(psver+1)
logit(1);

end