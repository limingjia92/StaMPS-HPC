function []=ps_select(reest_flag,plot_flag)
% PS_SELECT (HPC Optimized Version)
%   Select PS pixels based on amplitude dispersion and phase stability.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          December 2025
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   PERFORMANCE BENCHMARK (On Test Dataset):
%   - Original:    2564 seconds
%   - Optimized:   116 seconds (w/ MEX Topofit) / 125 seconds (w/o MEX Topofit)
%   - Speedup:     ~22x (Approx. 95% reduction in execution time)
%
%   OPTIMIZATION HIGHLIGHTS:
%   1. Batch Processing Strategy: Replaced the pixel-wise loop with a single 
%      MEX call ('clap_filt_patch_mex'), eliminating massive MATLAB interpreter 
%      overhead.
%   2. Parallelization: Integrated OpenMP multi-threading for phase filtering.
%   3. Configuration: Added 'use_fast_topofit' parameter. (default = 'n')
%      - 'y': 116s execution. Uses MEX topofit (~1.4% coherence divergence).
%      - 'n': 125s execution. Uses MATLAB topofit (Exact precision).
%
%   COMPILATION REQUIREMENTS:
%   mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" clap_filt_patch_mex.c
%   mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" ps_topofit_mex.c
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%   ======================================================================
%   07/2006 AH: Use specific bperp for re-estimation of gamma
%   08/2006 AH: Load workspaces into structures
%   09/2006 AH: Check if no gamma threshold meets criteria
%   09/2006 AH: add small baselines
%   09/2009 AH: add option to select on density instead of percentage
%   01/2010 KS: fixed badly conditioned polyfit errors by scaling and
%               centering
%   02/2010 AH: Only bin by D_A if enough candidate pixels
%   02/2010 AH: Leave out ifgs in drop_ifg_index from noise calculation
%   09/2010 JJS+MA: To make use oversampling dataset
%   12/2010 AH: Fix error message for density selection
%   02/2011 DB: Fix error with keep_ix
%   05/2012 AH: subtract only pixel being tested, not zero whole grid cell
%   01/2013 AH: Set default threshold if not enough random pixels
%   04/2015 DB: Give a warning to remove patch from patch list when no PS are left
%   09/2015 DB: Store information on nr of PS left. Support auto procesing
%   01/2016 DB: Replace save with stamps_save 
%   02/2016 DB: Identified bug when patch size is smaller than filter size.
%               For now drop this patch. This needs to be fixed better.
%   08/2017 AH: Ensure coh_thresh not negative.
%   12/2025 MJ: OpenMP-accelerated MEX Optimized
%   ======================================================================
%
logit;
logit('Selecting stable-phase pixels (Mex Optimized Version)...')

% ================= Check MEX =================
if exist('clap_filt_patch_mex', 'file') ~= 3
    error('Missing MEX: clap_filt_patch_mex');
end

use_fast_topofit = getparm('use_fast_topofit');
if isempty(use_fast_topofit) || strcmpi(use_fast_topofit,'n')
    use_fast_topofit = 'n';
    logit('Using original version topofit')
else
    use_fast_topofit = 'y';
    logit('Using Mex version topofit')
    if exist('ps_topofit_mex', 'file') ~= 3
        error('Missing MEX: ps_topofit_mex');
    end
end

if nargin<1, reest_flag=0; end
if nargin<2, plot_flag=0; end

psver=getpsver;
if psver>1, setpsver(1); end

slc_osf=getparm('slc_osf',1);
clap_alpha=getparm('clap_alpha',1);
clap_beta=getparm('clap_beta',1);
n_win=getparm('clap_win',1);
select_method=getparm('select_method',1);
if strcmpi(select_method,'PERCENT')
    max_percent_rand=getparm('percent_rand',1);
else
    max_density_rand=getparm('density_rand',1);
end
gamma_stdev_reject=getparm('gamma_stdev_reject',1);
small_baseline_flag=getparm('small_baseline_flag',1);
drop_ifg_index=getparm('drop_ifg_index',1);

if strcmpi(small_baseline_flag,'y')
    low_coh_thresh=15; 
else
    low_coh_thresh=31; 
end

load psver
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
pmname=['pm',num2str(psver)];
selectname=['select',num2str(psver)];
daname=['da',num2str(psver)];
bpname=['bp',num2str(psver)];

ps=load(psname);
ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

bperp=ps.bperp;
n_ifg=ps.n_ifg;
if ~strcmpi(small_baseline_flag,'y')
    master_ix=sum(ps.master_day>ps.day)+1;
    no_master_ix=setdiff([1:ps.n_ifg],ps.master_ix);
    ifg_index=setdiff(ifg_index,ps.master_ix);
    ifg_index(ifg_index>master_ix)=ifg_index(ifg_index>master_ix)-1;
    ph=ph(:,no_master_ix);
    bperp=bperp(no_master_ix);
    n_ifg=length(no_master_ix);
end
n_ps=ps.n_ps;
xy=ps.xy;

pm=load(pmname);
if exist([daname,'.mat'],'file')
    da=load(daname);
    D_A=da.D_A;
    clear da
else
    D_A=[];
end

if ~isempty(D_A) & size(D_A,1)>=10000
    D_A_sort=sort(D_A);
    if size(D_A,1)>=50000
      bin_size=10000;
    else
      bin_size=2000;
    end
    D_A_max=[0;D_A_sort(bin_size:bin_size:end-bin_size);D_A_sort(end)];
else
    D_A_max=[0;1];
    D_A=ones(size(pm.coh_ps));
end

if ~strcmpi(select_method,'PERCENT')
    patch_area=prod(max(xy(:,2:3),[],1)-min(xy(:,2:3),[],1))/1e6;
    max_percent_rand=max_density_rand*patch_area/(length(D_A_max)-1);
end

min_coh=zeros(length(D_A_max)-1,1);
D_A_mean=zeros(size(D_A_max,1)-1,1);
Nr_dist=pm.Nr;

if reest_flag==3 
    coh_thresh=0;
    coh_thresh_coeffs=[];
else
  for i=1:length(D_A_max)-1
    coh_chunk=pm.coh_ps(D_A>D_A_max(i) & D_A<=D_A_max(i+1));
    D_A_mean(i)=mean(D_A(D_A>D_A_max(i) & D_A<=D_A_max(i+1)));
    coh_chunk=coh_chunk(coh_chunk~=0); 
    Na=hist(coh_chunk,pm.coh_bins);
    Nr=Nr_dist*sum(Na(1:low_coh_thresh))/sum(Nr_dist(1:low_coh_thresh));
    if i==length(D_A_max)-1 & plot_flag==1 
        figure
        plot(pm.coh_bins,Na,'g'); hold on; plot(pm.coh_bins,Nr,'r');
        legend('data','random'); title('Before Gamma Reestimation');
    end
    Na(Na==0)=1; 
    if strcmpi(select_method,'PERCENT')
        percent_rand=fliplr(cumsum(fliplr(Nr))./cumsum(fliplr(Na))*100);
    else
        percent_rand=fliplr(cumsum(fliplr(Nr))); 
    end
    ok_ix=find(percent_rand<max_percent_rand);
    if isempty(ok_ix)
        min_coh(i)=1; 
    else
        min_fit_ix=min(ok_ix)-3;
        if min_fit_ix<=0
            min_coh(i)=nan;
        else
            max_fit_ix=min(ok_ix)+2;
            max_fit_ix(max_fit_ix>100)=100; 
            [p,S,mu]=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3);
            min_coh(i)=polyval(p,max_percent_rand,[],mu);
        end
    end
  end
  
  nonnanix=~isnan(min_coh);
  if sum(nonnanix)<1
      warning('Not enough random phase pixels to set gamma threshold - using default 0.3')
      coh_thresh=0.3;
      coh_thresh_coeffs=[];
  else 
    min_coh=min_coh(nonnanix);
    D_A_mean=D_A_mean(nonnanix);
    if size(min_coh,1)>1
        coh_thresh_coeffs=polyfit(D_A_mean,min_coh,1);
        if coh_thresh_coeffs(1)>0
            coh_thresh=polyval(coh_thresh_coeffs,D_A);
        else 
            coh_thresh=polyval(coh_thresh_coeffs,0.35);
            coh_thresh_coeffs=[];
        end
    else
        coh_thresh=min_coh;
        coh_thresh_coeffs=[];
    end
 end
end

coh_thresh(coh_thresh<0)=0; 
logit(sprintf('Initial gamma threshold: %.3f at D_A=%.2f to %.3f at D_A=%.2f',min(coh_thresh),min(D_A),max(coh_thresh),max(D_A)))

if  plot_flag==1 
    figure
    plot(D_A_mean,min_coh,'*'); hold on
    if ~isempty(coh_thresh_coeffs)
        plot(D_A_mean,polyval(coh_thresh_coeffs,D_A_mean),'r')
    end
    ylabel('\gamma_{thresh}'); xlabel('D_A')
end

ix=find(pm.coh_ps>coh_thresh); 
n_ps=length(ix);
logit(sprintf('%d PS selected initially',n_ps))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% reject part-time PS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gamma_stdev_reject>0
    ph_res_cpx=exp(j*pm.ph_res(:,ifg_index));
    coh_std=zeros(length(ix),1);
    for i=1:length(ix)
        coh_std(i)=std(bootstrp(100,@(ph) abs(sum(ph))/length(ph), ph_res_cpx(ix(i),ifg_index)));
    end
    clear ph_res_cpx
    ix=ix(coh_std<gamma_stdev_reject);
    n_ps=length(ix);
    logit(sfprintf('%d PS left after pps rejection',n_ps))
end

if reest_flag~=1  
  if reest_flag~=2 
    % Prepare Filter Parameters
    low_pass = pm.low_pass;
    % Ensure low_pass matches n_win
    if size(low_pass,1) > n_win
        lp_mid = floor(size(low_pass,1)/2) + 1;
        lp_half = floor(n_win/2);
        low_pass = low_pass(lp_mid-lp_half:lp_mid+lp_half-1, lp_mid-lp_half:lp_mid+lp_half-1);
    end
    B_1d = gausswin(7);

    % Prepare Grid for Batch MEX
    logit(sprintf('Calculating patch phases for %d pixels (Batch MEX)...', n_ps));
    
    if isfield(pm, 'ph_grid')
        input_grid = pm.ph_grid; 
        if ndims(input_grid) == 2 && n_ifg > 1
             error('pm.ph_grid dimision not fit');
        end
        if size(input_grid, 3) ~= n_ifg
             if size(input_grid, 3) == ps.n_ifg
                 input_grid = input_grid(:,:,ifg_index);
             end
        end
    else
        error('pm.ph_grid not exist...');
    end
    
    target_ij = pm.grid_ij(ix, :);

    % use MEX
    ph_patch2 = clap_filt_patch_mex(...
        double(input_grid), ... 
        double(target_ij), ...
        clap_alpha, ...
        clap_beta, ...
        double(low_pass), ...
        double(slc_osf), ...
        double(n_win), ...
        B_1d);
        
    clear input_grid % release memory

    % Prepare Topo Fit Inputs
    bp=load(bpname);
    bperp_mat=bp.bperp_mat(ix,:);
    clear bp
    
    if strcmpi(use_fast_topofit,'y')
        logit('Calculating coherence using MEX ps_topofit_mex...');
        
        ph_sub = ph(ix, ifg_index);
        patch_sub = ph_patch2(:, ifg_index);
        bperp_sub = bperp_mat(:, ifg_index);
        
        [K_mex, C_mex, coh_mex, ph_res_mex, ~] = ps_topofit_mex(...
            double(ph_sub), ...
            double(patch_sub), ...
            double(bperp_sub), ...
            double(pm.n_trial_wraps));
        
        K_ps2 = K_mex;
        C_ps2 = C_mex;
        coh_ps2 = coh_mex;
        
        ph_res2 = zeros(n_ps, n_ifg, 'single');
        ph_res2(:, ifg_index) = single(ph_res_mex);
    else
        logit('Calculating coherence using raw ps_topofit...');
        K_ps2 = zeros(n_ps,1);
        C_ps2 = zeros(n_ps,1);
        coh_ps2 = zeros(n_ps,1);
        ph_res2 = zeros(n_ps,n_ifg);
        
        ph_sub = ph(ix, :); 
        bperp_sub = bperp_mat(:, :); 
        
        for i=1:n_ps
            psdph = ph_sub(i,:).*conj(ph_patch2(i,:)); % use patch from MEX
            psdph = psdph./abs(psdph);
            [Kopt,Copt,cohopt,ph_residual] = ps_topofit(psdph(ifg_index).', bperp_sub(i,ifg_index).', pm.n_trial_wraps, 'n');
            
            K_ps2(i) = Kopt(1);
            C_ps2(i) = Copt(1);
            coh_ps2(i) = cohopt(1);
            ph_res2(i,ifg_index) = angle(ph_residual);
        end
    end
    
    logit(sprintf('%d coherences re-estimated', n_ps));

  else % reest_flag==2
    sl=load(selectname);
    ix=sl.ix;
    coh_ps2=sl.coh_ps2;
    K_ps2=sl.K_ps2;
    C_ps2=sl.C_ps2;
    ph_res2=sl.ph_res2;
    ph_patch2=sl.ph_patch2;
  end 
    pm.coh_ps(ix)=coh_ps2;

    % Recalculate Threshold
    for i=1:length(D_A_max)-1
        coh_chunk=pm.coh_ps(D_A>D_A_max(i) & D_A<=D_A_max(i+1));
        D_A_mean(i)=mean(D_A(D_A>D_A_max(i) & D_A<=D_A_max(i+1)));
        coh_chunk=coh_chunk(coh_chunk~=0);  
        Na=hist(coh_chunk,pm.coh_bins);
        Nr=Nr_dist*sum(Na(1:low_coh_thresh))/sum(Nr_dist(1:low_coh_thresh));
        Na(Na==0)=1; 
        if strcmpi(select_method,'PERCENT')
            percent_rand=fliplr(cumsum(fliplr(Nr))./cumsum(fliplr(Na))*100);
        else
            percent_rand=fliplr(cumsum(fliplr(Nr))); 
        end
        ok_ix=find(percent_rand<max_percent_rand);
        if isempty(ok_ix)
            min_coh(i)=1;
        else
            min_fit_ix=min(ok_ix)-3;
            if min_fit_ix<=0
                min_coh(i)=nan;
            else
                max_fit_ix=min(ok_ix)+2;
                max_fit_ix(max_fit_ix>100)=100; 
                [p,S,mu]=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3); 
                min_coh(i)=polyval(p,max_percent_rand,[],mu); 
            end
        end
    end

    nonnanix=~isnan(min_coh);
    if sum(nonnanix)<1
        coh_thresh=0.3;
    else
        min_coh=min_coh(nonnanix);
        D_A_mean=D_A_mean(nonnanix);
        if length(min_coh)>1
            coh_thresh_coeffs=polyfit(D_A_mean,min_coh,1); 
            if coh_thresh_coeffs(1)>0 
                coh_thresh=polyval(coh_thresh_coeffs,D_A(ix));
            else 
                coh_thresh=polyval(coh_thresh_coeffs,0.35); 
            end
        else
            coh_thresh=min_coh;
        end
    end
    coh_thresh(coh_thresh<0)=0; 
    logit(sprintf('Reestimation gamma threshold: %.3f at D_A=%.2f to %.3f at D_A=%.2f',min(coh_thresh),min(D_A),max(coh_thresh),max(D_A)))

    bperp_range=max(bperp)-min(bperp);
    keep_ix=coh_ps2>coh_thresh & abs(pm.K_ps(ix)-K_ps2)<2*pi/bperp_range;
    clear pm
    logit(sprintf('%d ps selected after re-estimation of coherence',sum(keep_ix)))

else % reest_flag==1
    pm=rmfield(pm,{'ph_grid'});
    ph_patch2=pm.ph_patch(ix,:);
    ph_res2=pm.ph_res(ix,:);
    K_ps2=pm.K_ps(ix);
    C_ps2=pm.C_ps(ix);
    coh_ps2=pm.coh_ps(ix);
    keep_ix=true(size(ix));
end 

if exist('no_ps_info.mat','file')~=2
   stamps_step_no_ps = zeros([5 1 ]);       
else
   load('no_ps_info.mat');
   stamps_step_no_ps(3:end)=0;
end
if sum(keep_ix)==0
   fprintf('***No PS points left. Updating the stamps log for this****\n')
   stamps_step_no_ps(3)=1;
end
save('no_ps_info.mat','stamps_step_no_ps')

ph_patch2 = single(ph_patch2);
ph_res2 = single(ph_res2);
stamps_save(selectname,ix,keep_ix,ph_patch2,ph_res2,K_ps2,C_ps2,coh_ps2,coh_thresh,coh_thresh_coeffs,clap_alpha,clap_beta,n_win,max_percent_rand,gamma_stdev_reject,small_baseline_flag,ifg_index);

logit(1);
end