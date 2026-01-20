function logit(logmsg, whereto)
%LOGIT (HPC Optimized Version)
%   Write message to STAMPS.log and/or stdout
%   logit(logmsg, whereto)
%
%   INPUTS:
%       logmsg  - Message string OR 0 ('Starting') OR 1 ('Finished')
%       whereto - 0 (Default): Write to STAMPS.log AND stdout
%                 1: Write to STAMPS.log ONLY
%                 2: Write to stdout ONLY
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          January 2026
%   Version:       1.0
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, March 2010
%

    if nargin < 1
        logmsg = 0;
    end
    if nargin < 2
        whereto = 0;
    end

    if isnumeric(logmsg)
       switch logmsg
           case 0
               logmsg = 'Starting';
           case 1
               logmsg = 'Finished';
           otherwise
               logmsg = num2str(logmsg); % Fallback for other numbers
       end
    end

    if length(logmsg) > 1 && strcmp(logmsg(end-1:end), '\n')
        logmsg = logmsg(1:end-2);
    end

    [ST] = dbstack(1); 
    if ~isempty(ST)
        fname = upper(ST(1).name);
    else
        fname = 'COMMAND_LINE';
    end

    t_now = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'); 
    t_str = char(t_now); % Convert to char for efficient fprintf
        
    % CASE: Write to File (whereto 0 or 1)
    if whereto == 0 || whereto == 1
        logid = fopen('STAMPS.log', 'a');
        if logid > 0
            fprintf(logid, '%s %-16s %s\n', t_str, fname, logmsg);
            fclose(logid);
        end
    end

    % CASE: Write to Screen (whereto 0 or 2)
    if whereto == 0 || whereto == 2
        fprintf('%s: %s\n', fname, logmsg);
    end

end