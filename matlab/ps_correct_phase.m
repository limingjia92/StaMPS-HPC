function []=ps_correct_phase
% PS_CORRECT_PHASE() correct phase from estimate of look angle error
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   07/2006 AH: Use specific bperp for correction
%   09/2006 AH: add small baselines 
%   12/2025 Mingjia: bug fix and code optimized
%   ==========================================================
logit;
fprintf('Correcting phase for look angle error...\n')

small_baseline_flag=getparm('small_baseline_flag',1);

load psver
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
pmname=['pm',num2str(psver)];
rcname=['rc',num2str(psver)];
bpname=['bp',num2str(psver)];

ps=load(psname);
pm=load(pmname);
bp=load(bpname);

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

% Ensure parameters are single column vectors
K_ps=single(pm.K_ps);
C_ps=single(pm.C_ps);
master_ix=sum(ps.master_day>ps.day)+1;

if strcmpi(small_baseline_flag,'y')
    % OPTIMIZATION: Implicit expansion
    phase_correction = K_ps .* single(bp.bperp_mat); % subtract range error 
    ph_rc = ph .* exp(-1j * phase_correction); 
    
    stamps_save(rcname,ph_rc);
else    
    % Construct bperp_mat with single precision zeros
    % Ensure input parts are single to avoid double creation
    bperp_part1 = single(bp.bperp_mat(:,1:ps.master_ix-1));
    bperp_part2 = single(bp.bperp_mat(:,ps.master_ix:end));
    
    % Concatenate
    bperp_mat = [bperp_part1, zeros(ps.n_ps,1,'single'), bperp_part2];
    
    % OPTIMIZATION: Implicit expansion
    phase_term = K_ps .* bperp_mat + C_ps;
    ph_rc = ph .* exp(-1j * phase_term);
    
    % Prepare re-referenced phase
    ph_patch_single = single(pm.ph_patch);
    ph_reref = [ph_patch_single(:,1:master_ix-1), ones(ps.n_ps,1,'single'), ph_patch_single(:,master_ix:end)];
    
    stamps_save(rcname,ph_rc,ph_reref);
end

logit(1);

end