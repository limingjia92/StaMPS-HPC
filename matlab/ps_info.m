function []=ps_info(force_single_master)
%PS_INFO Display information on the dataset (PS and SBAS).
%   ps_info(force_single_master)
%
%   Usage:
%       ps_info()  - Auto-detects mode (SBAS pairs or PS single master).
%       ps_info(1) - Forces single-master display mode.
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
%   1. Unified Interface: Merged legacy sb_info logic into ps_info.
%   2. Dependency Bypass: Removed forced calls to ps_calc_ifg_std.
%   3. Output Formatting: Modernized dates to ISO 8601 (YYYY-MM-DD) and 
%      added dynamic scaling for pixel counts (e.g., millions) to 
%      facilitate automated HPC log parsing and quick human readability.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, March 2007
%   ======================================================================

% 1. Initialization and parameter handling
if nargin < 1
    force_single_master = 0;
end

load psver
psname=['ps',num2str(psver)];
ifgstdname=['ifgstd',num2str(psver)];

small_baseline_flag=getparm('small_baseline_flag');
is_sbas = strcmpi(small_baseline_flag,'y');
ps=load(psname);

% 2. Load standard deviation passively
has_ifg_std = false;
if exist([ifgstdname,'.mat'],'file')
    stdin=load(ifgstdname);
    ifg_std=stdin.ifg_std;
    has_ifg_std = true;
end

if is_sbas && force_single_master == 0
    % 3. Display SBAS Mode (Interferogram Pairs)
    if isfield(ps,'ifgday')
        fprintf('\n--------------------------------------------------------------------------\n');
        fprintf(' ID   Reference Date   to  Secondary Date      B_perp    Noise (ifg_std)\n');
        fprintf('--------------------------------------------------------------------------\n');
        
        for i=1:ps.n_ifg
            date_m = datestr(ps.ifgday(i,1), 'yyyy-mm-dd');
            date_s = datestr(ps.ifgday(i,2), 'yyyy-mm-dd');
            bperp_str = sprintf('%5d m', round(ps.bperp(i)));
            
            if has_ifg_std && ifg_std(i) ~= 0
                fprintf('%3d   %s   to   %s   %s   %8.3f deg\n', i, date_m, date_s, bperp_str, ifg_std(i));
            else
                fprintf('%3d   %s   to   %s   %s\n', i, date_m, date_s, bperp_str);
            end
        end
        fprintf('--------------------------------------------------------------------------\n');
    else
        fprintf('\n[Warning] Data flagged as SBAS, but no "ifgday" field found in ps.mat.\n');
        return;
    end
    
else
    % 4. Display PS Mode (Single Master Sequence)
    if is_sbas
        G=zeros(ps.n_ifg,ps.n_image);
        for i=1:ps.n_ifg
             G(i,ps.ifgday_ix(i,1))=-1;
             G(i,ps.ifgday_ix(i,2))=1;
        end
        G=G(:,[1:ps.master_ix-1,ps.master_ix+1:end]);
        bperp=[G\double(ps.bperp)];
        bperp=[bperp(1:ps.master_ix-1);0;bperp(ps.master_ix:end)];
    else
        bperp=ps.bperp;
    end

    fprintf('\n------------------------------------------------------\n');
    fprintf(' ID       Date         B_perp       Noise (ifg_std)\n');
    fprintf('------------------------------------------------------\n');

    for i=1:size(ps.day,1)
        date_str = datestr(ps.day(i), 'yyyy-mm-dd');
        bperp_str = sprintf('%5d m', round(bperp(i)));
        
        if has_ifg_std && ifg_std(i) ~= 0
            fprintf('%3d   %s   %s   %8.3f deg\n', i, date_str, bperp_str, ifg_std(i));
        else
            fprintf('%3d   %s   %s\n', i, date_str, bperp_str);
        end 
    end
    fprintf('------------------------------------------------------\n');
end

% 5. Print Shared Footer (PS Pixel Count)
n_ps = ps.n_ps;
if n_ps >= 1e6
    ps_str = sprintf('%.2f million', n_ps / 1e6);
elseif n_ps >= 1e3
    ps_str = sprintf('%.1f thousand', n_ps / 1e3);
else
    ps_str = sprintf('%d', n_ps);
end

exact_str = regexprep(num2str(n_ps), '(?<=\d)(?=(\d{3})+(?!\d))', ',');

fprintf('Number of stable-phase pixels: %s (%s)\n', ps_str, exact_str);
if is_sbas && force_single_master == 0
    fprintf('--------------------------------------------------------------------------\n\n');
else
    fprintf('------------------------------------------------------\n\n');
end