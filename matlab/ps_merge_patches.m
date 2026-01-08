function []=ps_merge_patches(psver)
%PS_MERGE_PATCHES (HPC Optimized Version)
%   Merge overlapping patches into a single dataset using Parallel I/O.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          December 2025
%   Version:       2.0 (Parallel-Variable-Stream)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   PERFORMANCE BENCHMARK (On Test Dataset):
%   - Wall Clock:  Reduced from ~29 min to ~10.5 min (~2.8x Speedup).
%   - CPU Efficiency: Total CPU cycles reduced by ~91% (11.6x Efficiency Gain).
%     (Reduced User Time from 16,516s to 1,423s).
%   - Memory Usage: Stable at ~15.8GB (Identical to original).
%     (Achieved via Variable-by-Variable processing strategy).
%
%   OPTIMIZATION HIGHLIGHTS:
%   1. Two-Phase Architecture: 
%      - Phase 1: Lightweight serial scan to calculate indices & global sort order.
%      - Phase 2: Parallel processing of heavy data variables.
%   2. Variable-Centric Parfor: Instead of looping patches (serial), we loop 
%      variables. Inside each variable task, patches are processed in parallel (`parfor`).
%   3. Cell Array Buffering: Replaced dynamic matrix concatenation with Cell Arrays 
%      to solve Parfor slicing issues and eliminate memory fragmentation.
%   4. Selective I/O: Only loads the specific variable being processed, drastically 
%      reducing I/O overhead compared to loading entire workspaces.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, September 2006
%   ======================================================================
%   01/2007 AH: scla and scn added
%   01/2007 AH: patch.list added
%   06/2009 AH: file existence checks to search current directory only
%   06/2009 AH: move mean amplitude merge to end to save memory
%   06/2009 AH: only save scla and scn when present to save memory
%   09/2009 AH: add option to resample to coarser sampling 
%   09/2009 AH: reduce memory needs further
%   11/2009 AH: ensure mean amplitude width is always correct
%   06/2010 AH: estimate weights directly from residuals
%   10/2010 DB: Fix when patch_noover does not have PS (resampling)
%   10/2010 DB: Fix when PS all rejected when sum weight<min_weight (resampling)
%   06/2010 AH: Move mean amplitude merging to ps_load_mean_amp.m 
%   02/2011 DB: Fix dimension for min computation of ps.xy 
%   09/2015 AH: Delete previous merged amplitude files
%   06/2017 DB: Include stamps save for large variables
%   10/2017 DB: If inc angle file is present also merge it.
%   ======================================================================

logit;
fprintf('Merging patches ...\n')

%% 1. Initial Setup
if nargin < 1
    psver=2;
end
patch_list_file = 'patch.list';

if exist(['./',patch_list_file],'file')
    dirname=struct;
    fid=fopen(patch_list_file,'r');
    i=0;
    while feof(fid)==0
        i=i+1;
        dirname(i).name=fgetl(fid);
    end
    fclose(fid);
else
    dirname=dir('PATCH_*');
end

if ~exist('./parms.mat','file')
    disp('Use the parms.mat in the 1st PATCH dir')
    copyfile(strcat(dirname(1).name,filesep,'parms.mat'),'.');
end

pwdpath=pwd;
small_baseline_flag=getparm('small_baseline_flag');
grid_size=getparm('merge_resample_size',1);
merge_stdev=getparm('merge_standard_dev',1);
phase_accuracy=10*pi/180;
min_weight=1/merge_stdev^2;
randn('state',1001);
max_coh=abs(sum(exp(1i*randn(1000,1)*phase_accuracy)))/1000;

% Variable Names
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];
bpname=['bp',num2str(psver)];
laname=['la',num2str(psver)];
incname=['inc',num2str(psver)];
hgtname=['hgt',num2str(psver)];
scnname=['scn',num2str(psver)];
sclaname=['scla',num2str(psver)];
sclasbname=['scla_sb',num2str(psver)];
phuwname=['phuw',num2str(psver)];

n_patch=length(dirname);

%% 2. PHASE 1: Index Calculation (Serial - Light IO)
fprintf('--- Phase 1: Calculating Merge Indices ---\n');

% Pre-allocate MetaInfo
MetaInfo(n_patch) = struct('ix', [], 'n_ps_g', 0, 'ps_weight', [], 'ps_snr', [], 'f_ix', [], 'l_ix', []);

ij_global = zeros(0,2); 
lonlat_global = zeros(0,2); 
coh_ps_global = zeros(0,1,'single'); 
remove_ix_global = logical([]);

for k=1:n_patch
    if ~isempty(dirname(k).name)
        fprintf('   Analyzing structure of %s\n', dirname(k).name);
        cd(dirname(k).name);
        
        ps = load(psname);
        n_ifg = ps.n_ifg; 
        
        patch_ij = load('patch_noover.in');
        ix = (ps.ij(:,2)>=patch_ij(3)-1 & ps.ij(:,2)<=patch_ij(4)-1 & ...
              ps.ij(:,3)>=patch_ij(1)-1 & ps.ij(:,3)<=patch_ij(2)-1);
        
        if sum(ix)==0, ix_no_ps = 1; else, ix_no_ps = 0; end
        
        if grid_size == 0
            [~,IA,IB] = intersect(ps.ij(ix,2:3), ij_global, 'rows');
            remove_ix_global = [remove_ix_global; IB]; 
            
            [~,IA,~] = intersect(ps.ij(:,2:3), ij_global, 'rows');
            ix_ex = true(ps.n_ps, 1);
            ix_ex(IA) = 0;
            ix(ix_ex) = 1;
            
            MetaInfo(k).ix = ix;
            MetaInfo(k).n_ps_g = sum(ix); 
            
            ij_global = [ij_global; ps.ij(ix,2:3)];
            lonlat_global = [lonlat_global; ps.lonlat(ix,:)];
            
            % Minimal load for coh_ps
            if exist([pmname, '.mat'], 'file')
                pm_tmp = load(pmname, 'coh_ps');
                if isfield(pm_tmp, 'coh_ps')
                     coh_ps_global = [coh_ps_global; pm_tmp.coh_ps(ix)];
                else
                     coh_ps_global = [coh_ps_global; zeros(sum(ix),1,'single')];
                end
            else
                coh_ps_global = [coh_ps_global; zeros(sum(ix),1,'single')];
            end
            
        elseif grid_size ~= 0 && ix_no_ps ~= 1
            % Resampling Logic
            clear g_ij
            xy_min = min(ps.xy(ix,:), 1); 
            g_ij(:,1) = ceil((ps.xy(ix,3)-xy_min(3)+1e-9)/grid_size);
            g_ij(:,2) = ceil((ps.xy(ix,2)-xy_min(2)+1e-9)/grid_size);
            
            [g_ij, ~, g_ix] = unique(g_ij, 'rows');
            [g_ix, sort_ix_local] = sort(g_ix);
            
            ix_indices = find(ix);
            ix_indices = ix_indices(sort_ix_local); 
            
            pm = load(pmname, 'ph_res', 'coh_ps', 'C_ps');
            pm.ph_res=angle(exp(1i*(pm.ph_res-repmat(pm.C_ps,1,size(pm.ph_res,2)))));
            if small_baseline_flag~='y'
                pm.ph_res=[pm.ph_res,pm.C_ps];
            end
            
            sigsq_noise=var([pm.ph_res],0,2); 
            coh_ps_all=abs(sum(exp(1i*[pm.ph_res]),2))/n_ifg;
            coh_ps_all(coh_ps_all>max_coh)=max_coh; 
            sigsq_noise(sigsq_noise<phase_accuracy^2)=phase_accuracy^2; 
            
            ps_weight = single(1./sigsq_noise(ix_indices)); 
            ps_snr = single(1./(1./coh_ps_all(ix_indices).^2 - 1));
            clear pm
            
            l_ix_local = [find(diff(g_ix)); size(g_ix,1)];
            f_ix_local = [1; l_ix_local(1:end-1)+1];
            n_ps_g = size(f_ix_local, 1);
            
            weightsave = zeros(n_ps_g, 1, 'single'); 
            mask_ix = true(length(ix_indices), 1);
            
            for i=1:n_ps_g
                weights = ps_weight(f_ix_local(i):l_ix_local(i));
                weightsave(i) = sum(weights);
                if weightsave(i) < min_weight
                    mask_ix(f_ix_local(i):l_ix_local(i)) = false;
                end
            end
            
            valid_groups = weightsave >= min_weight;
            
            MetaInfo(k).ix = ix_indices(mask_ix); 
            MetaInfo(k).n_ps_g = sum(valid_groups);
            MetaInfo(k).ps_weight = ps_weight(mask_ix);
            MetaInfo(k).ps_snr = ps_snr(mask_ix);
            
            g_ix = g_ix(mask_ix);
            l_ix_new = [find(diff(g_ix)); size(g_ix,1)];
            f_ix_new = [1; l_ix_new(1:end-1)+1];
            
            MetaInfo(k).f_ix = f_ix_new;
            MetaInfo(k).l_ix = l_ix_new;
            
            ij_g = zeros(MetaInfo(k).n_ps_g, 2); 
            lonlat_g = zeros(MetaInfo(k).n_ps_g, 2); 
            coh_ps_g = zeros(MetaInfo(k).n_ps_g, 1, 'single');
            
            ps_ij_sel = ps.ij(MetaInfo(k).ix, :);       
            ps_lonlat_sel = ps.lonlat(MetaInfo(k).ix, :);
            
            for i=1:MetaInfo(k).n_ps_g
                 w = repmat(MetaInfo(k).ps_weight(f_ix_new(i):l_ix_new(i)), 1, 2);
                 sub_ij = ps_ij_sel(f_ix_new(i):l_ix_new(i), 2:3);
                 sub_ll = ps_lonlat_sel(f_ix_new(i):l_ix_new(i), :);
                 
                 ij_g(i,:) = round(sum(sub_ij.*w, 1) ./ sum(w(:,1))); 
                 lonlat_g(i,:) = sum(sub_ll.*w, 1) ./ sum(w(:,1));
                 
                 w_coh = MetaInfo(k).ps_weight(f_ix_new(i):l_ix_new(i));
                 snr_val = sqrt(sum(w_coh.^2, 1));
                 coh_ps_g(i) = sqrt(1./(1+1./snr_val));
            end
            
            ij_global = [ij_global; ij_g];
            lonlat_global = [lonlat_global; lonlat_g];
            coh_ps_global = [coh_ps_global; coh_ps_g];
        end
        cd(pwdpath);
    end
end

n_ps_total_raw = size(ij_global, 1);
fprintf('Total points (raw): %d\n', n_ps_total_raw);

%% 3. Global Sort & Clean
fprintf('--- Calculating Global Coordinates & Removing Duplicates ---\n');

keep_ix = true(n_ps_total_raw, 1);
keep_ix(remove_ix_global) = 0;
coh_ps_weed = coh_ps_global(keep_ix);

lonlat_filtered = lonlat_global(keep_ix,:); 

[dummy, I] = unique(lonlat_filtered, 'rows');
dups = setxor(I, [1:size(lonlat_filtered,1)]'); 
keep_ix_num = find(keep_ix); 

if ~isempty(dups)
    fprintf('   Resolving %d duplicate pixel groups...\n', length(dups));
    for i=1:length(dups)
        target_lat = lonlat_filtered(dups(i),1);
        target_lon = lonlat_filtered(dups(i),2);
        
        dups_ix_weed = find(lonlat_filtered(:,1)==target_lat & lonlat_filtered(:,2)==target_lon);
        dups_ix = keep_ix_num(dups_ix_weed);
        
        [~, max_I] = max(coh_ps_weed(dups_ix_weed));
        keep_ix(dups_ix([1:end]~=max_I)) = 0; 
    end
end

final_indices = find(keep_ix);
final_lonlat = lonlat_global(final_indices, :); 

ll0 = (max(final_lonlat) + min(final_lonlat)) / 2;
xy = llh2local(final_lonlat', ll0) * 1000; 
xy = xy';

heading = getparm('heading');
if isempty(heading), heading=0; end
theta = (180-heading)*pi/180;
if theta>pi, theta=theta-2*pi; end
rotm = [cos(theta), sin(theta); -sin(theta), cos(theta)];
xynew = (rotm * xy')';

if (max(xynew(:,1))-min(xynew(:,1)) < max(xy(:,1))-min(xy(:,1))) && ...
   (max(xynew(:,2))-min(xynew(:,2)) < max(xy(:,2))-min(xy(:,2)))
   xy = xynew;
   disp(['   Rotating xy by ', num2str(theta*180/pi), ' degrees']);
end

xy = single(xy);
[~, sort_order] = sortrows(xy, [2,1]);
xy = xy(sort_order, :);
xy = round(xy * 1000) / 1000; 

final_indices_sorted = final_indices(sort_order);
n_ps = length(final_indices_sorted);

fprintf('   Writing ps%d.mat (contains %d pixels)\n', psver, n_ps);

% Save PS struct
ps_template_path = [dirname(1).name, filesep, psname];
if exist([ps_template_path, '.mat'], 'file')
    ps_new = load(ps_template_path);
else
    ps_new = struct(); 
    ps_new.n_ifg = n_ifg;
end

ps_new.n_ps = n_ps;
ps_new.ij = [[1:n_ps]', ij_global(final_indices_sorted,:)];
ps_new.xy = [[1:n_ps]', xy]; 
ps_new.lonlat = final_lonlat(sort_order, :);
ps_new.ll0 = ll0;

save(psname, '-struct', 'ps_new');
clear ps_new ij_global lonlat_global coh_ps_global keep_ix lonlat_filtered xy

%% 4. PHASE 2: Variable Processing (Parallel & Memory Optimized)
Tasks = {
    'bp',  0, bpname;
    'la',  0, laname;
    'inc', 0, incname;
    'hgt', 0, hgtname;
    'ph',  0, phname; 
    'rc',  1, rcname;
    'pm',  0, pmname;
    'phuw',0, phuwname;
    'scla',1, sclaname;
    'scla_sb',1, sclasbname;
    'scn', 0, scnname;
};

% User specified list for Double Conversion
vars_to_double = {'hgt', 'inc', 'la', 'ph', 'pm', 'rc'};

% Get dimensions
cd(dirname(1).name);
if exist([bpname, '.mat'], 'file'), dummy_bp = load(bpname); n_cols_bp = size(dummy_bp.bperp_mat, 2); else n_cols_bp=0; end
if exist([phname, '.mat'], 'file'), dummy_ph = load(phname); n_cols_ifg = size(dummy_ph.ph, 2); else n_cols_ifg=0; end
cd(pwdpath);

% Start Parallel Pool
pool = gcp('nocreate');
if isempty(pool)
    parpool; 
end

PatchNames = {dirname.name};

for t = 1:size(Tasks, 1)
    varType = Tasks{t, 1};
    saveName = Tasks{t, 3};
    
    % Check existence in first patch
    fileExists = exist([dirname(1).name, filesep, saveName, '.mat'], 'file');
    if ~fileExists && ~strcmp(varType, 'inc') && ~strcmp(varType, 'rc')
        continue;
    end

    if fileExists
        fprintf('--- Processing Variable: %s ---\n', saveName);
    else
        fprintf('--- Creating Empty Variable: %s ---\n', saveName);
    end

    % =========================================================================
    % CASE 1: PM (Complex Structure - 5 Sub-variables)
    % =========================================================================
    if strcmp(varType, 'pm') && fileExists
        
        % Initialize CELL ARRAYS for parfor output
        OUT_PH_PATCH = cell(1, n_patch);
        OUT_PH_RES   = cell(1, n_patch);
        OUT_K_PS     = cell(1, n_patch);
        OUT_C_PS     = cell(1, n_patch);
        OUT_COH_PS   = cell(1, n_patch);
        
        parfor k = 1:n_patch
            if MetaInfo(k).n_ps_g == 0, continue; end
            
            loaded = load([PatchNames{k}, filesep, saveName, '.mat']);
            info = MetaInfo(k);

            % Process Logic
            raw_patch = loaded.ph_patch(info.ix, :);
            if isfield(loaded,'ph_res'), raw_res=loaded.ph_res(info.ix,:); else, raw_res=raw_patch*0; end
            if isfield(loaded,'K_ps'), raw_K=loaded.K_ps(info.ix,:); else, raw_K=zeros(sum(info.ix),1); end
            if isfield(loaded,'C_ps'), raw_C=loaded.C_ps(info.ix,:); else, raw_C=zeros(sum(info.ix),1); end
            
            if grid_size == 0
                OUT_PH_PATCH{k} = raw_patch;
                OUT_PH_RES{k}   = raw_res;
                OUT_K_PS{k}     = raw_K;
                OUT_C_PS{k}     = raw_C;
                OUT_COH_PS{k}   = loaded.coh_ps(info.ix, :);
            else
                res_patch = zeros(info.n_ps_g, n_cols_ifg, 'single');
                res_res   = zeros(info.n_ps_g, n_cols_ifg, 'single');
                res_K     = zeros(info.n_ps_g, 1, 'single');
                res_C     = zeros(info.n_ps_g, 1, 'single');
                res_coh   = zeros(info.n_ps_g, 1, 'single');
                
                for i=1:info.n_ps_g
                    w_snr = repmat(info.ps_snr(info.f_ix(i):info.l_ix(i)), 1, n_cols_ifg);
                    res_patch(i,:) = sum(raw_patch(info.f_ix(i):info.l_ix(i),:) .* w_snr, 1);
                    res_res(i,:)   = sum(raw_res(info.f_ix(i):info.l_ix(i),:)   .* w_snr, 1);
                    
                    snr_sq_sum = sum(w_snr(:,1).^2, 1);
                    res_coh(i) = sqrt(1 ./ (1 + 1 ./ sqrt(snr_sq_sum)));
                    
                    w_var = info.ps_weight(info.f_ix(i):info.l_ix(i));
                    sum_w_var = sum(w_var);
                    res_K(i) = sum(raw_K(info.f_ix(i):info.l_ix(i)) .* w_var) ./ sum_w_var;
                    res_C(i) = sum(raw_C(info.f_ix(i):info.l_ix(i)) .* w_var) ./ sum_w_var;
                end
                OUT_PH_PATCH{k} = res_patch;
                OUT_PH_RES{k}   = res_res;
                OUT_K_PS{k}     = res_K;
                OUT_C_PS{k}     = res_C;
                OUT_COH_PS{k}   = res_coh;
            end
            loaded = [];
        end

        % Vertcat (Fast), Slice, Convert, Save
        % Note: vertcat of cells preserves order 1..n_patch, aligning with Phase 1.
        
        ph_patch = vertcat(OUT_PH_PATCH{:}); clear OUT_PH_PATCH;
        ph_patch = double(ph_patch(final_indices_sorted, :));
        
        ph_res   = vertcat(OUT_PH_RES{:}); clear OUT_PH_RES;
        ph_res   = double(ph_res(final_indices_sorted, :));
        
        K_ps     = vertcat(OUT_K_PS{:}); clear OUT_K_PS;
        K_ps     = double(K_ps(final_indices_sorted, :));
        
        C_ps     = vertcat(OUT_C_PS{:}); clear OUT_C_PS;
        C_ps     = double(C_ps(final_indices_sorted, :));
        
        coh_ps   = vertcat(OUT_COH_PS{:}); clear OUT_COH_PS;
        coh_ps   = double(coh_ps(final_indices_sorted, :));
        
        stamps_save(saveName, ph_patch, ph_res, K_ps, C_ps, coh_ps);
        clear ph_patch ph_res K_ps C_ps coh_ps

    % =========================================================================
    % CASE 2: RC / SCLA (Multi-variable)
    % =========================================================================
    elseif strcmp(varType, 'rc') || strcmp(varType, 'scla') || strcmp(varType, 'scla_sb')
        
        is_rc = strcmp(varType, 'rc');
        OUT_MAIN = cell(1, n_patch);
        OUT_SEC  = cell(1, n_patch);
        OUT_TRD  = cell(1, n_patch);
        
        parfor k = 1:n_patch
             if MetaInfo(k).n_ps_g == 0, continue; end
             info = MetaInfo(k);
             
             if fileExists
                loaded = load([PatchNames{k}, filesep, saveName, '.mat']);
             else
                loaded = struct(); 
             end

             % -- Handle Main Var --
             if fileExists
                 if is_rc, f='ph_rc'; else, f='ph_scla'; end
                 raw_d = loaded.(f)(info.ix, :);
                 
                 if grid_size == 0
                     OUT_MAIN{k} = raw_d;
                 else
                     res_d = zeros(info.n_ps_g, n_cols_ifg, 'single');
                     for i=1:info.n_ps_g
                        if is_rc
                             w = repmat(info.ps_snr(info.f_ix(i):info.l_ix(i)), 1, n_cols_ifg);
                             res_d(i,:) = sum(raw_d(info.f_ix(i):info.l_ix(i),:) .* w, 1);
                        else
                             w = repmat(info.ps_weight(info.f_ix(i):info.l_ix(i)), 1, n_cols_ifg);
                             res_d(i,:) = sum(raw_d(info.f_ix(i):info.l_ix(i),:) .* w, 1) ./ sum(w(:,1));
                        end
                     end
                     OUT_MAIN{k} = res_d;
                 end
             elseif is_rc
                 % RC special case: if not exists, create zeros
                 OUT_MAIN{k} = zeros(info.n_ps_g, n_cols_ifg, 'single');
             end
             
             % -- Handle Sec Var --
             if is_rc 
                 % RC: ph_reref
                 if ~strcmpi(small_baseline_flag, 'y') || (~fileExists && ~strcmpi(small_baseline_flag, 'y'))
                     if fileExists && isfield(loaded, 'ph_reref')
                         raw_d = loaded.ph_reref(info.ix, :);
                     else
                         raw_d = zeros(sum(info.ix), n_cols_ifg, 'single');
                     end
                     
                     if grid_size==0
                         OUT_SEC{k} = raw_d;
                     else
                         res_d = zeros(info.n_ps_g, n_cols_ifg, 'single');
                         for i=1:info.n_ps_g
                              w = repmat(info.ps_snr(info.f_ix(i):info.l_ix(i)), 1, n_cols_ifg);
                              res_d(i,:) = sum(raw_d(info.f_ix(i):info.l_ix(i),:) .* w, 1);
                         end
                         OUT_SEC{k} = res_d;
                     end
                 end
             else 
                 % SCLA: K_ps_uw, C_ps_uw
                 if fileExists
                     raw_k = loaded.K_ps_uw(info.ix, :);
                     raw_c = loaded.C_ps_uw(info.ix, :);
                     if grid_size==0
                         OUT_SEC{k} = raw_k;
                         OUT_TRD{k} = raw_c;
                     else
                         res_k = zeros(info.n_ps_g, 1, 'single');
                         res_c = zeros(info.n_ps_g, 1, 'single');
                         for i=1:info.n_ps_g
                             w = info.ps_weight(info.f_ix(i):info.l_ix(i));
                             res_k(i) = sum(raw_k(info.f_ix(i):info.l_ix(i)).*w) ./ sum(w);
                             res_c(i) = sum(raw_c(info.f_ix(i):info.l_ix(i)).*w) ./ sum(w);
                         end
                         OUT_SEC{k} = res_k;
                         OUT_TRD{k} = res_c;
                     end
                 end
             end
        end
        
        if is_rc
            VAR_MAIN = vertcat(OUT_MAIN{:}); clear OUT_MAIN;
            if ~isempty(VAR_MAIN)
                 ph_rc = double(VAR_MAIN(final_indices_sorted, :)); 
                 mask = ph_rc~=0; ph_rc(mask) = ph_rc(mask) ./ abs(ph_rc(mask));
            else
                 ph_rc = [];
            end
            clear VAR_MAIN
            
            VAR_SEC = vertcat(OUT_SEC{:}); clear OUT_SEC;
            if ~isempty(VAR_SEC)
                ph_reref = double(VAR_SEC(final_indices_sorted, :));
            else
                ph_reref = [];
            end
            clear VAR_SEC
            stamps_save(saveName, ph_rc, ph_reref);
            clear ph_rc ph_reref
        else
            ph_scla = vertcat(OUT_MAIN{:}); clear OUT_MAIN;
            ph_scla = single(ph_scla(final_indices_sorted, :));
            
            K_ps_uw = vertcat(OUT_SEC{:}); clear OUT_SEC;
            K_ps_uw = single(K_ps_uw(final_indices_sorted, :));
            
            C_ps_uw = vertcat(OUT_TRD{:}); clear OUT_TRD;
            C_ps_uw = single(C_ps_uw(final_indices_sorted, :));
            
            stamps_save(saveName, ph_scla, K_ps_uw, C_ps_uw);
            clear ph_scla K_ps_uw C_ps_uw
        end

    % =========================================================================
    % CASE 3: BP (BPERP)
    % =========================================================================
    elseif strcmp(varType, 'bp') && fileExists
         OUT_BP = cell(1, n_patch);
         
         parfor k=1:n_patch
             if MetaInfo(k).n_ps_g == 0, continue; end
             info = MetaInfo(k);
             
             loaded = load([PatchNames{k}, filesep, saveName, '.mat']);
             raw_bp = loaded.bperp_mat(info.ix, :);
             
             if grid_size==0
                 OUT_BP{k} = raw_bp;
             else
                 res_bp = zeros(info.n_ps_g, n_cols_bp, 'single');
                 for i=1:info.n_ps_g
                      w = repmat(info.ps_weight(info.f_ix(i):info.l_ix(i)), 1, n_cols_bp);
                      w(w==0) = 1e-9;
                      res_bp(i,:) = sum(raw_bp(info.f_ix(i):info.l_ix(i),:) .* w, 1) ./ sum(w(:,1));
                 end
                 OUT_BP{k} = res_bp;
             end
         end
         
         bperp_mat = vertcat(OUT_BP{:}); clear OUT_BP;
         bperp_mat = single(bperp_mat(final_indices_sorted, :)); 
         stamps_save(saveName, bperp_mat);
         clear bperp_mat

    % =========================================================================
    % CASE 4: Simple Variables (ph, inc, la, hgt, scn, phuw)
    % =========================================================================
    elseif fileExists || strcmp(varType, 'inc') 
        
        % Identify field name
        if strcmp(varType, 'ph'), f='ph'; n_cols=n_cols_ifg;
        elseif strcmp(varType, 'phuw'), f='ph_uw'; n_cols=n_cols_ifg;
        elseif strcmp(varType, 'scn'), f='ph_scn_slave'; n_cols=n_cols_ifg;
        elseif strcmp(varType, 'la') || strcmp(varType, 'inc') || strcmp(varType, 'hgt'), f=varType; n_cols=1;
        end
        
        OUT_VAR = cell(1, n_patch);
        
        % Use a flag to avoid checking fileExists inside parfor
        do_process = false;
        if fileExists, do_process = true; end

        parfor k=1:n_patch
            if MetaInfo(k).n_ps_g == 0, continue; end
            
            % Special handle for INC which might not exist on disk but needs 0s
            if ~do_process && strcmp(varType, 'inc')
                 OUT_VAR{k} = zeros(MetaInfo(k).n_ps_g, 1, 'single');
                 continue;
            end
            
            info = MetaInfo(k);
            loaded = load([PatchNames{k}, filesep, saveName, '.mat']);
            raw_d = loaded.(f)(info.ix, :);
            
            if grid_size == 0
                OUT_VAR{k} = raw_d;
            else
                res_d = zeros(info.n_ps_g, n_cols, 'single');
                for i=1:info.n_ps_g
                     if strcmp(varType, 'ph')
                         w = repmat(info.ps_snr(info.f_ix(i):info.l_ix(i)), 1, n_cols);
                         res_d(i,:) = sum(raw_d(info.f_ix(i):info.l_ix(i),:) .* w, 1);
                     else
                         w = repmat(info.ps_weight(info.f_ix(i):info.l_ix(i)), 1, n_cols);
                         res_d(i,:) = sum(raw_d(info.f_ix(i):info.l_ix(i),:) .* w, 1) ./ sum(w(:,1));
                     end
                end
                OUT_VAR{k} = res_d;
            end
        end
        
        SINGLE_VAR = vertcat(OUT_VAR{:}); clear OUT_VAR;
        
        if isempty(SINGLE_VAR)
            val = [];
        else
            if ismember(varType, vars_to_double)
                val = double(SINGLE_VAR(final_indices_sorted, :));
            else
                val = SINGLE_VAR(final_indices_sorted, :);
            end
        end
        clear SINGLE_VAR;
        
        eval([f ' = val;']);
        args = {f};
        cmd_str = 'stamps_save(saveName'; 
        for i = 1:length(args)
            cmd_str = [cmd_str, ', ', args{i}];
        end
        cmd_str = [cmd_str, ');'];
        eval(cmd_str);
        
        eval(['clear ', f]);
    end
end

save psver psver
logit(1);
end