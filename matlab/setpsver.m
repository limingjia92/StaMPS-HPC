function setpsver(new_psver)
% SETPSVER Set the processing version (step) and save to psver.mat
%
%   Andy Hooper, June 2006
%   Revised by Mingjia, January 2026

    if nargin < 1 || ~isnumeric(new_psver)
        error('Error: setpsver requires a numeric input argument.');
    end
    
    ver_file = 'psver.mat';
    
    current_ver = getpsver(); 
    fprintf('psver currently: %d\n', current_ver);

    psver = new_psver; 
    
    try
        save(ver_file, 'psver');
        fprintf('psver now set to:  %d\n', psver);
    catch ME
        error('Failed to write to %s. Reason: %s', ver_file, ME.message);
    end
end