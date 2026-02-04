function []=setparm(parmname,value,newflag)
%SETPARM (HPC Optimized Version)
%   sets parameter value in parms.mat 
%   USAGE:
%       setparm(parmname, value)          % Update an existing parameter
%       setparm(parmname, value, 1)       % Create a NEW parameter
%       setparm(parmname, [], -1)         % DELETE a parameter
%       setparm(parmname, nan)            % Reset parameter to default
%       setparm()                         % Display all parameters
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0 (HPC-Math)
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, June 2006
%

    parmfile_name = 'parms.mat';
    parmfile_path = parmfile_name; % Default to current dir

    % --- 1. Load Parameters (Priority: Current Dir > Parent Dir) ---
    if exist(parmfile_name,'file')
        parms = load(parmfile_name);
    elseif exist(['..' filesep parmfile_name],'file')
        parmfile_path = ['..' filesep parmfile_name];
        parms = load(parmfile_path);
    else
        parms = struct(); % Create new if nothing exists
    end

    % --- 2. Input Handling ---
    if nargin == 0
        % Display parameters
        disp('--- Current Parameters ---');
        disp(['Source: ' parmfile_path]);
        disp(orderfields(parms));
        return;
    end

    if nargin < 2
        error('Format is: SETPARM(PARMNAME,VALUE,[NEWFLAG])');
    end

    % --- 3. Parameter Name Resolution ---
    % If NOT creating a new parameter (newflag~=1), we must find the full name
    if nargin <= 2 || newflag ~= 1
        if ~isfield(parms, parmname)
            valid_fields = fieldnames(parms);
            match_idx = find(strncmpi(parmname, valid_fields, length(parmname)));
            
            if length(match_idx) > 1
                error(['Parameter ''',parmname,'*'' is not unique. Matches: ', strjoin(valid_fields(match_idx)', ', ')])
            elseif isempty(match_idx)
                if nargin > 2 && newflag == -1
                     fprintf('Parameter %s does not exist, nothing to delete.\n', parmname);
                     return;
                else
                     error(['Parameter ''',parmname,''' does not exist. Use newflag=1 to create it.'])
                end
            end
            parmname = valid_fields{match_idx};
        end
    end

    % --- 4. Logging Helper ---
    % Helper string for log
    val_str = '[]';
    if isnumeric(value)
        val_str = num2str(value);
    elseif ischar(value)
        val_str = value;
    end

    % --- 5. Action Execution ---
    
    % CASE A: Delete Parameter (newflag == -1)
    if nargin > 2 && newflag == -1
        if isfield(parms, parmname)
            parms = rmfield(parms, parmname);
            save(parmfile_path, '-struct', 'parms');
            logit([parmname, ' removed from ', parmfile_path], 0);
        end
        return;
    end

    % CASE B: Reset to Default (value is NaN)
    if isnumeric(value) && isscalar(value) && isnan(value)
        if strcmpi(parmname,'small_baseline_flag')
            error('Default reset option not possible for small_baselines_flag');
        end
        
        % Remove current value
        if isfield(parms, parmname)
            parms = rmfield(parms, parmname);
            save(parmfile_path, '-struct', 'parms'); % Save state before loading defaults
        end
        
        % Reload defaults (this script injects defaults into parms.mat if missing)
        ps_parms_default; 
        
        % Read back the restored value for display
        new_val = getparm(parmname); 
        disp([parmname, ' reset to default value: ', num2str(new_val)]);
        return;
    end

    % CASE C: Set/Update Parameter (Standard)
    % Log the change
    logit([parmname, ' = ', val_str], 0);
    
    % Update struct
    parms.(parmname) = value;
    
    % Save to disk
    save(parmfile_path, '-struct', 'parms');

end