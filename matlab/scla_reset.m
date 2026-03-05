function [] = scla_reset(patches_flag)
%SCLA_RESET Reset estimated SCLA and master atmosphere/orbit error files.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization & Fixes:
%   1. Resource Management: Added missing 'fclose(fid)' to prevent file 
%      handle leaks during massive batch processing on clusters.
%   2. Code Refactoring: Replaced repetitive if-delete blocks with a clean 
%      cell array loop for better readability and maintainability.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, Jan 2007
%   ======================================================================

logit;

if nargin < 1
    patches_flag = 'y';
end

logit('Resetting SCLA error and master atmosphere/orbit error...', 2)

i = 0;
patchdir = struct('name', {});

% --- Find patch directories ---
if strcmpi(patches_flag, 'y')
    if exist('patch.list', 'file')
        fid = fopen('patch.list');
        while 1
            nextline = fgetl(fid);
            if ischar(nextline)
                i = i + 1;
                patchdir(i).name = nextline;
            else
                break;
            end
        end
        fclose(fid); % [BUG FIX] 
    else
        patchdir = dir('PATCH_*');
        i = length(patchdir);
    end
end

% Include current directory
patchdir(i+1).name = '.';
currdir = pwd;

% --- Delete SCLA files in all target directories ---
for j = 1:length(patchdir)
    if ~isempty(patchdir(j).name)
        cd(patchdir(j).name)
        logit([pwd, ' reset'])

        if exist('./psver.mat', 'file')
            load psver
            
            % Target files to clean up
            files_to_delete = {
                ['scla', num2str(psver), '.mat'], ...
                ['scla_smooth', num2str(psver), '.mat'], ...
                ['scla_sb', num2str(psver), '.mat'], ...
                ['scla_smooth_sb', num2str(psver), '.mat']
            };
            
            % Refactored deletion loop
            for k = 1:length(files_to_delete)
                if exist(files_to_delete{k}, 'file')
                    delete(files_to_delete{k});
                end
            end
        end
        cd(currdir)
    end
end

logit(1);
end