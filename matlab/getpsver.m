function psver = getpsver()
% GETPSVER Get the current processing version (step) from psver.mat
%
%   Andy Hooper, June 2006
%   Revised by Mingjia, January 2026

    ver_file = 'psver.mat';

    if exist(ver_file, 'file')
        data = load(ver_file, 'psver');
        if isfield(data, 'psver')
            psver = data.psver;
        else
            warning('File %s exists but does not contain "psver" variable. Defaulting to 0.', ver_file);
            psver = 1;
        end
    else
        psver = 1; 
    end
end