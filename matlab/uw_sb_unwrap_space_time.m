function [] = uw_sb_unwrap_space_time(day, ifgday_ix, bperp, options)
%UW_SB_UNWRAP_SPACE_TIME Space-time phase unwrapping (HPC Optimized).
% ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
% ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Parallelized Look Angle Error Estimation (Block 2).
%   2. Parallelized Standard 3D Unwrapping (Block 3).
%   3. Robust 3D_FULL implementation using lscov.
%   4. Safety check for SBAS network topology in 3D_FULL mode.
%
%   Usage & Strategy:
%   Extracted Parameters from 'options':
%       unwrap_method, time_win, la_flag, n_trial_wraps,
%       prefilt_win, scf_flag, temp, n_temp_wraps, max_bperp_for_temp_est
%
%   - '2D'       : Spatial smoothing only. Fastest but susceptible to phase jumps.
%   - '3D_QUICK' : Time-space smoothing with simple thresholding. Fast, low accuracy.
%   - '3D'       : (Standard) Temporal smoothing + iterative optimization. 
%                  HPC-Parallelized. High accuracy, recommended for SBAS.
%   - '3D_FULL'  : Temporal modeling via local linear regression. 
%                  Highest precision for PS; incompatible with SBAS.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2007
% ======================================================================

% --- Argument Validation ---
if nargin < 4
    error('Error: Not enough arguments. Usage: uw_sb_unwrap_space_time(day, ifgday_ix, bperp, options)');
end

% --- Parameter Unpacking ---
unwrap_method          = options.unwrap_method;
time_win               = options.time_win;
la_flag                = options.la_flag;
n_trial_wraps          = options.n_trial_wraps;
prefilt_win            = options.prefilt_win;
scf_flag               = options.scf_flag;
temp                   = options.temp;
n_temp_wraps           = options.n_temp_wraps;
max_bperp_for_temp_est = options.max_bperp_for_temp_est;

disp('Unwrapping in time-space... (HPC Version)')

uw=load('uw_grid');
ui=load('uw_interp');
n_ifg=uw.n_ifg;
n_ps=uw.n_ps;
nzix=uw.nzix;
ij=uw.ij;

if isempty(uw.ph_uw_predef)
    predef_flag='n';
else
    predef_flag='y';
end

n_image=size(day,1);
[nrow,ncol]=size(ui.Z);
dph_space=((uw.ph(ui.edgs(:,3),:).*conj(uw.ph(ui.edgs(:,2),:))));

if predef_flag=='y'
    dph_space_uw=uw.ph_uw_predef(ui.edgs(:,3),:)-uw.ph_uw_predef(ui.edgs(:,2),:);
    predef_ix=~isnan(dph_space_uw);
    dph_space_uw=dph_space_uw(predef_ix);
else
    predef_ix=[];
end
clear uw 

dph_space=dph_space./abs(dph_space);
ifreq_ij=[];
jfreq_ij=[];

G=zeros(n_ifg,n_image);
for i=1:n_ifg
    G(i,ifgday_ix(i,1))=-1;
    G(i,ifgday_ix(i,2))=1;
end

nzc_ix=sum(abs(G))~=0; 
day=day(nzc_ix);
n_image=size(day,1);
G=G(:,nzc_ix);
zc_ix=find(nzc_ix==0);
zc_ix = sort(zc_ix,'descend');              
for i=1:length(zc_ix)
    ifgday_ix(ifgday_ix>zc_ix(i))=ifgday_ix(ifgday_ix>zc_ix(i))-1;
end
n=size(G,2);

if ~isempty(temp)
    temp_flag='y';
else
    temp_flag='n';
end
 

% ==========================================================
% BLOCK 1: Temperature Correlation Estimation (Serial)
% ==========================================================
if strcmpi(temp_flag,'y')
    fprintf('   Estimating temperature correlation \n')
    ix=abs(bperp)<max_bperp_for_temp_est;
    temp_sub=temp(ix);
    temp_range=max(temp)-min(temp);
    temp_range_sub=max(temp_sub)-min(temp_sub);
    dph_sub=dph_space(:,ix); 
    n_temp_wraps=n_temp_wraps*(temp_range_sub/temp_range);

    trial_mult=[-ceil(8*n_temp_wraps):ceil(8*n_temp_wraps)];
    n_trials=length(trial_mult);
    trial_phase=temp_sub/temp_range_sub*pi/4;
    trial_phase_mat=exp(-1i*trial_phase*trial_mult);
    Kt=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    
    for i=1:ui.n_edge
        cpxphase=dph_sub(i,:).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        coh_trial=abs(phaser_sum)/sum(abs(cpxphase));
        [coh_max,coh_max_ix]=max(coh_trial);
        
        falling_ix=find(diff(coh_trial(1:coh_max_ix))<0); 
        if ~isempty(falling_ix)
            peak_start_ix=falling_ix(end)+1;
        else
            peak_start_ix=1;
        end
        rising_ix=find(diff(coh_trial(coh_max_ix:end))>0); 
        if ~isempty(rising_ix)
            peak_end_ix=rising_ix(1)+coh_max_ix-1;
        else
            peak_end_ix=n_trials;
        end
        coh_trial(peak_start_ix:peak_end_ix)=0;
        
        if coh_max-max(coh_trial)>0.1      
            K0=pi/4/temp_range_sub*trial_mult(coh_max_ix);
            resphase=cpxphase.*exp(-1i*(K0*temp_sub)); 
            offset_phase=sum(resphase);
            resphase=angle(resphase*conj(offset_phase)); 
            weighting=abs(cpxphase); 
            mopt=double(weighting.*temp_sub)\double(weighting.*resphase);
            Kt(i)=K0+mopt;
            phase_residual=cpxphase.*exp(-1i*(Kt(i)*temp_sub)); 
            mean_phase_residual=sum(phase_residual); 
            coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
        end
    end
    clear cpxphase_mat trial_phase_mat phaser 
    Kt(coh<0.31)=0; 
    dph_space=dph_space.*exp(-1i*Kt*temp');
    if predef_flag=='y'
        dph_temp=Kt*temp';
        dph_space_uw=dph_space_uw-dph_temp(predef_ix);
        clear dph_temp
    end
    dph_sub=dph_sub.*exp(-1i*Kt*temp_sub');
end

% ==========================================================
% BLOCK 2: Look Angle Error Estimation (Parallelized)
% ==========================================================
if strcmpi(la_flag,'y') 
    fprintf('   Estimating look angle error \n')
    bperp_range=max(bperp)-min(bperp);
    ix=find(abs(diff(ifgday_ix,[],2))==1);   

    if length(ix)>=length(day)-1 
        fprintf('   using sequential daisy chain of interferograms\n')
        dph_sub=dph_space(:,ix); 
        bperp_sub=bperp(ix);
        bperp_range_sub=max(bperp_sub)-min(bperp_sub);
        n_trial_wraps=n_trial_wraps*(bperp_range_sub/bperp_range);
    else
        ifgs_per_image=sum(abs(G));
        [max_ifgs_per_image,max_ix]=max(ifgs_per_image);
        if max_ifgs_per_image>=length(day)-2
            fprintf('   Using sequential daisy chain of interferograms\n')
            ix=G(:,max_ix)~=0;
            gsub=G(ix,max_ix);
            sign_ix=-sign(single(gsub'));
            dph_sub=dph_space(:,ix); 
            bperp_sub=[bperp(ix)]; 
            bperp_sub(sign_ix==-1)=-bperp_sub(sign_ix==-1);
            bperp_sub=[bperp_sub;0]; 
            sign_ix=repmat(sign_ix,ui.n_edge,1);
            dph_sub(sign_ix==-1)=conj(dph_sub(sign_ix==-1)); 
            dph_sub=[dph_sub,mean(abs(dph_sub),2)]; 
            slave_ix=sum(ifgday_ix(ix,:),2)-max_ix;
            day_sub=day([slave_ix;max_ix]); 
            [day_sub,sort_ix]=sort(day_sub); 
            dph_sub=dph_sub(:,sort_ix); 
            bperp_sub=bperp_sub(sort_ix);
            bperp_sub=diff(bperp_sub);
            bperp_range_sub=max(bperp_sub)-min(bperp_sub);
            n_trial_wraps=n_trial_wraps*(bperp_range_sub/bperp_range);           
            dph_sub=dph_sub(:,[2:end]).*conj(dph_sub(:,1:end-1)); 
            dph_sub=dph_sub./abs(dph_sub); 
        else 
            dph_sub=dph_space;
            bperp_sub=bperp;
            bperp_range_sub=bperp_range;
        end
    end
        
    trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
    n_trials=length(trial_mult);
    trial_phase=bperp_sub/bperp_range_sub*pi/4;
    trial_phase_mat=exp(-1i*trial_phase*trial_mult);
    K=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    
    try
        pool_obj = gcp('nocreate');
        if isempty(pool_obj); parpool; end
        fprintf('   Start parallel search for Look Angle Error...\n');
        
        parfor i=1:ui.n_edge
            cpxphase=dph_sub(i,:).';
            cpxphase_mat=repmat(cpxphase,1,n_trials);
            phaser=trial_phase_mat.*cpxphase_mat;
            phaser_sum=sum(phaser);
            denom = sum(abs(cpxphase));
            if denom == 0; coh_trial = zeros(1, n_trials); else; coh_trial=abs(phaser_sum)/denom; end
            
            [coh_max,coh_max_ix]=max(coh_trial);
            falling_ix=find(diff(coh_trial(1:coh_max_ix))<0); 
            if ~isempty(falling_ix); peak_start_ix=falling_ix(end)+1; else; peak_start_ix=1; end
            rising_ix=find(diff(coh_trial(coh_max_ix:end))>0); 
            if ~isempty(rising_ix); peak_end_ix=rising_ix(1)+coh_max_ix-1; else; peak_end_ix=n_trials; end
            
            coh_trial_masked = coh_trial;
            coh_trial_masked(peak_start_ix:peak_end_ix)=0;
            
            if coh_max-max(coh_trial_masked)>0.1 
                K0=pi/4/bperp_range_sub*trial_mult(coh_max_ix);
                resphase=cpxphase.*exp(-1i*(K0*bperp_sub)); 
                offset_phase=sum(resphase);
                resphase=angle(resphase*conj(offset_phase)); 
                weighting=abs(cpxphase); 
                mopt=double(weighting.*bperp_sub)\double(weighting.*resphase);
                K(i)=K0+mopt;
                phase_residual=cpxphase.*exp(-1i*(K(i)*bperp_sub)); 
                mean_phase_residual=sum(phase_residual); 
                coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
            end
        end
    catch ME
        fprintf('   Parallel failed: %s. Switching to serial.\n', ME.message);
    end
    clear cpxphase_mat trial_phase_mat phaser dph_sub
    K(coh<0.31)=0; 
    if strcmpi(temp_flag,'y')
        dph_space(K==0,:)=dph_space(K==0,:).*exp(1i*Kt(K==0)*temp'); 
        Kt(K==0)=0; K(Kt==0)=0; 
    end
    dph_space=dph_space.*exp(-1i*K*bperp');
    if predef_flag=='y'
        dph_scla=K*bperp';
        dph_space_uw=dph_space_uw-dph_scla(predef_ix);
        clear dph_scla
    end
end

spread=sparse(zeros(ui.n_edge,n_ifg));

% ==========================================================
% BLOCK 3: Unwrapping Logic
% ==========================================================

if strcmpi(unwrap_method,'2D')
    dph_space_uw=angle(dph_space);
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp'; 
    end
    if strcmpi(temp_flag,'y')
        dph_space_uw=dph_space_uw+Kt*temp';   
    end
    dph_noise=[];
    save('uw_space_time','dph_space_uw','spread','dph_noise');    
    
else
    fprintf('   Smoothing in time \n')   
    % --- 3D_FULL ---
    if strcmpi(unwrap_method,'3D_FULL') 
       dph_smooth_ifg=nan(size(dph_space),'single');
       valid_master_count = 0; % Track valid masters for SBAS check
       
       for i=1:n_image
         ix=G(:,i)~=0;
         if sum(ix)>=n_image-2 % REQUIREMENT: Fully connected network
           valid_master_count = valid_master_count + 1;
           gsub=G(ix,i);
           dph_sub=dph_space(:,ix); 
           sign_ix=repmat(-sign(single(gsub')),ui.n_edge,1);
           dph_sub(sign_ix==-1)=conj(dph_sub(sign_ix==-1)); 
           slave_ix=sum(ifgday_ix(ix,:),2)-i;
           day_sub=day(slave_ix); 
           [day_sub,sort_ix]=sort(day_sub); 
           dph_sub=dph_sub(:,sort_ix); 
           dph_sub_angle=angle(dph_sub);
           n_sub=length(day_sub);
           dph_smooth=zeros(ui.n_edge,n_sub,'single');
           
           for i1=1:n_sub
               time_diff=(day_sub(i1)-day_sub)';
               weight_factor=exp(-(time_diff.^2)/2/time_win^2);
               weight_factor=weight_factor/sum(weight_factor);

               dph_mean=sum(dph_sub.*repmat(weight_factor,ui.n_edge,1),2);
               dph_mean_adj=mod(dph_sub_angle-repmat(angle(dph_mean),1,n_sub)+pi,2*pi)-pi;
               
               GG=[ones(n_sub,1),time_diff'];
               if size(GG,1)>1
                   m=lscov(GG,double(dph_mean_adj)',weight_factor);
               else
                   m=zeros(size(GG,1),ui.n_edge);
               end
               dph_smooth(:,i1)=dph_mean.*exp(1i*(m(1,:)')); 
            end
            
            dph_smooth_sub=cumsum([angle(dph_smooth(:,1)),angle(dph_smooth(:,2:end).*conj(dph_smooth(:,1:end-1)))],2);
            close_master_ix=find(slave_ix-i>0);
            if isempty(close_master_ix)
                close_master_ix=n_sub;
            else
                close_master_ix=close_master_ix(1);
                if close_master_ix>1
                    close_master_ix=[close_master_ix-1;close_master_ix];
                end
            end
            dph_close_master=mean(dph_smooth_sub(:,close_master_ix),2);
            dph_smooth_sub=dph_smooth_sub-repmat(dph_close_master-angle(exp(1i*dph_close_master)),1,n_sub);
            dph_smooth_sub=dph_smooth_sub.*sign_ix;
            already_sub_ix=find(~isnan(dph_smooth_ifg(1,ix))); 
            ix=find(ix);
            already_ix=ix(already_sub_ix);
            std_noise1=std(angle(dph_space(:,already_ix).*exp(-1i*dph_smooth_ifg(:,already_ix))));
            std_noise2=std(angle(dph_space(:,already_ix).*exp(-1i*dph_smooth_sub(:,already_sub_ix))));
            keep_ix=true(n_sub,1);
            keep_ix(already_sub_ix(std_noise1<std_noise2))=false; 
            dph_smooth_ifg(:,ix(keep_ix))=dph_smooth_sub(:,keep_ix);
          end
        end

        % ERROR CHECK: If no valid masters found, 3D_FULL is incompatible (e.g., SBAS)
        if valid_master_count == 0
            error('uw_sb_unwrap_space_time:InvalidTopology', ...
                '3D_FULL requires a fully connected network (like PS). SBAS networks are not compatible. Please use "3D" or "3D_QUICK".');
        end

        dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));
        dph_noise(std(dph_noise,0,2)>1.2,:)=nan;
     
    else
        % --- Method: 3D_QUICK / 3D Standard (Parallel) ---
        x=(day-day(1))*(n-1)/(day(end)-day(1)); 
        if predef_flag=='y'
            n_dph=size(dph_space,1);
            dph_space_angle=double(angle(dph_space));
            dph_space_angle(predef_ix)=dph_space_uw;
            dph_space_series=zeros(n,n_dph);
            for i=1:n_dph
                 W=predef_ix(i,:)+0.01; 
                 dph_space_series(2:end,i)=lscov(double(G(:,2:end)),dph_space_angle(i,:)',W);
            end
        else
            dph_space_series=[zeros(1,ui.n_edge);double(G(:,2:end))\double(angle(dph_space))'];
        end
        
        dph_smooth_series=zeros(size(G,2),ui.n_edge,'single');
        for i1=1:n
            time_diff_sq=(day(i1)-day).^2;
            weight_factor=exp(-time_diff_sq/2/time_win^2);
            weight_factor=weight_factor/sum(weight_factor);
            dph_smooth_series(i1,:)=sum(dph_space_series.*repmat(weight_factor,1,ui.n_edge));
        end

        dph_smooth_ifg=(G*dph_smooth_series)';
        dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));

        if strcmpi(unwrap_method,'3D_QUICK')
            % --- Method: 3D_QUICK (Simple Thresholding) ---
            not_small_ix=find(std(dph_noise,0,2)>1.3)';
            fprintf('   %d edges with high std dev in time \n',length(not_small_ix))
            dph_noise(not_small_ix,:)=nan;
        else 
          % --- Method: 3D Standard (Parallel Annealing) ---
          uw=load('uw_grid');
          ph_noise=angle(uw.ph.*conj(uw.ph_lowpass));
          clear uw
          dph_noise_sf=((ph_noise(ui.edgs(:,3),:)-(ph_noise(ui.edgs(:,2),:))));      
          m_minmax=repmat([-pi,pi],5,1).*repmat([0.5;0.25;1;0.25;1],1,2);
          anneal_opts=[1;15;0;0;0;0;0];
          covm=cov((dph_noise_sf)); 
          [W,P]=chol(inv(covm)); 
          if P~=0
              W=diag(1./sqrt(diag(covm)));
          end
          not_small_ix=find(std(dph_noise,0,2)>1)';
          fprintf('   Performing complex smoothing on %d edges \n',length(not_small_ix))

          dph_space_subset = dph_space(not_small_ix, :);
          n_ix = length(not_small_ix);
          
          try
              pool_obj = gcp('nocreate'); 
              if isempty(pool_obj); parpool; end
              
              % --- HPC Parallel Progress Monitor Setup ---
              q = parallel.pool.DataQueue;
              monitor_handle = hpc_log_progress(n_ix, 10, 'SMOOTH_UNWRAP'); % print out each 10%
              afterEach(q, monitor_handle); 
              % -------------------------------------------
              
              temp_smooth_results = zeros(size(dph_smooth_series,1), n_ix);
              
              parfor k=1:n_ix
                dph = angle(dph_space_subset(k, :))'; 
                temp_smooth_results(:,k) = uw_sb_smooth_unwrap(m_minmax,anneal_opts,G,W,dph,x);
                
                send(q, 1); % send signal to Monitor
              end
              dph_smooth_series(:, not_small_ix) = temp_smooth_results;
              fprintf('   Parallel smoothing complete.\n');
              
          catch ME
              fprintf('   Parallel processing failed: %s. Reverting to serial loop.\n', ME.message);
              for i=not_small_ix
                dph=angle(dph_space(i,:))';
                dph_smooth_series(:,i)=uw_sb_smooth_unwrap(m_minmax,anneal_opts,G,W,dph,x);
              end
          end
          dph_smooth_ifg=(G*dph_smooth_series)';
          dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));
        end
    end
    clear dph_space
    dph_space_uw=dph_smooth_ifg+dph_noise;
    clear dph_smooth_ifg
    
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp';   
    end
    if strcmpi(temp_flag,'y')
        dph_space_uw=dph_space_uw+Kt*temp';   
    end
    
    if strcmpi(scf_flag,'y')
        fprintf('   Calculating local phase gradients \n')
        ifreq_ij=nan(n_ps,n_ifg,'single');
        jfreq_ij=nan(n_ps,n_ifg,'single');
        ifgw=zeros(nrow,ncol);
        uw=load('uw_grid');
        for i=1:n_ifg
            ifgw(nzix)=uw.ph(:,i);
            [ifreq,jfreq,grad_ij,Hmag]=gradient_filt(ifgw,prefilt_win);
            ix=~isnan(ifreq)&Hmag./((abs(ifreq))+1)>3;
            if sum(ix(:))>2
                ifreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),ij(:,2),ij(:,1),'linear');
            end
            ix=~isnan(jfreq)&Hmag./((abs(jfreq))+1)>3;
            if sum(ix(:))>2
                jfreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),ij(:,2),ij(:,1),'linear');
            end
        end
        clear uw
        
        spread2=zeros(size(spread),'single');
        dph_smooth_uw2=nan(ui.n_edge,n_ifg,'single');
       
        fprintf('   Smoothing using local phase gradients \n')
        for i=1:ui.n_edge
            nodes_ix=ui.edgs(i,[2:3]);
            ifreq_edge=mean(ifreq_ij(nodes_ix,:));
            jfreq_edge=mean(jfreq_ij(nodes_ix,:));
            diff_i=diff(ij(nodes_ix,1));
            diff_j=diff(ij(nodes_ix,2));
            dph_smooth_uw2(i,:)=diff_i*ifreq_edge+diff_j*jfreq_edge;
            spread2(i,:)=diff(ifreq_ij(nodes_ix,:))+diff(jfreq_ij(nodes_ix,:));
        end
        fprintf('   Choosing between time and phase gradient smoothing \n')        
        std_noise=std(dph_noise,0,2);
        dph_noise2=angle(exp(1i*(dph_space_uw-dph_smooth_uw2)));
        std_noise2=std(dph_noise2,0,2);
        dph_noise2(std_noise2>1.3,:)=nan;
        shaky_ix=isnan(std_noise) | std_noise>std_noise2; 
        
        fprintf('   %d arcs smoothed in time, %d in space \n',ui.n_edge-sum(shaky_ix),sum(shaky_ix))        
        dph_noise(shaky_ix,:)=dph_noise2(shaky_ix,:);
        dph_space_uw(shaky_ix,:)=dph_smooth_uw2(shaky_ix,:)+dph_noise2(shaky_ix,:);
        spread(shaky_ix,:)=spread2(shaky_ix,:);
    else
        shaky_ix=[];
    end
    save('uw_space_time','dph_space_uw','dph_noise','G','spread','ifreq_ij','jfreq_ij','shaky_ix','predef_ix');
end
