function [value, parmname] = getparm(parmname, printflag)
%GETPARM (HPC Optimized Version)
%   gets parameter value from parms.mat 
%   USAGE:
%       [value, full_name] = getparm(parmname)
%       getparm()  % Displays all parameters
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, July 2006
%

    % --- 1. Load Parameters (Priority: Current Dir > Parent Dir) ---
    parmfile_name = 'parms.mat';
    
    if exist(['.' filesep parmfile_name], 'file')
        parms = load(parmfile_name);
    elseif exist(['..' filesep parmfile_name], 'file')
        parms = load(['..' filesep parmfile_name]);
    else
        % Requirement 5: Only run defaults if file is missing
        fprintf('StaMPS-HPC: parms.mat not found. Initializing defaults...\n');
        ps_parms_default(); 
        
        % Try loading again after creation
        if exist(parmfile_name, 'file')
            parms = load(parmfile_name);
        else
            error('StaMPS-HPC: Critical Error. Could not create or load parms.mat.');
        end
    end

    % --- 2. Viewer Mode (Requirement 3) ---
    if nargin < 1
        disp('--- Current Parameters ---');
        disp(orderfields(parms));
        return;
    end

    if nargin < 2
        printflag = 0;
    end

    % --- 3. Parameter Name Resolution (Requirement 2: Fuzzy Matching) ---
    % Check if exact match exists first to save time
    if ~isfield(parms, parmname)
        valid_fields = fieldnames(parms);
        % Case-insensitive, length-based partial match
        match_idx = find(strncmpi(parmname, valid_fields, length(parmname)));
        
        if length(match_idx) > 1
            % Formatting error message for ambiguous matches
            matches = strjoin(valid_fields(match_idx)', ', ');
            error(['Parameter ''%s*'' is not unique. Matches: %s'], parmname, matches);
        elseif isempty(match_idx)
            % Parameter genuinely doesn't exist
            value = [];
            parmname = [];
            fprintf('Parameter ''%s'' not found in parms.mat. Returning empty.\n', parmname);
            return;
        else
            % Unique match found, update parmname to full name
            parmname = valid_fields{match_idx};
        end
    end

    % --- 4. Retrieve Value (Requirement 6: Dynamic Field Access) ---
    if isfield(parms, parmname)
        value = parms.(parmname);
    else
        value = [];
    end

    % --- 5. Logging / Display ---
    if printflag ~= 0
        if isnumeric(value)
            % Format numeric array (handle scalar or vector)
            val_str = num2str(value(:)'); 
            if length(val_str) > 200
                val_str = [val_str(1:200) '...']; % Truncate long arrays
            end
            msg = sprintf('%s = %s', parmname, val_str);
        elseif ischar(value)
            msg = sprintf('%s = ''%s''', parmname, value);
        else
            msg = sprintf('%s = [Complex Structure]', parmname);
        end
        logit(msg);
    end
end