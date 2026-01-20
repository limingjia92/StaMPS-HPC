function [] = stamps_save(save_name, varargin)
% STAMPS_SAVE (HPC Compatible Version)
%  save resuts for StaMPS
% INPUTS: 
% save_name     String with the filename of the save datafile
% varagin       Variables which need to be save. 
%               Note these are the actual variabels and not their names.
%
% By Bekaert David - Jan 2016
% modifications:
%
% DB    12/2016     Update to the version of TRAIN which allows for
%                   appending, error catching and folder generation if
%                   needed. Refer to TRAIN toolbox for the latest version.
% AH    08/2017     If no directory specified, ensure only current directory 
%                   searched. Add '.mat' if not specified.
% MJ    01/2026     Enforces '-v7.3' (HDF5) for all NEW files. 
%                   And Uses structures instead of 'eval' for safety and speed.

if ~strcmpi(save_name(max(1,end-3):end), '.mat')
    save_name = [save_name, '.mat'];
end

[file_path, ~, ~] = fileparts(save_name);
if ~isempty(file_path)
    if ~exist(file_path, 'dir')
        mkdir(file_path);
    end
else
    save_name = ['./', save_name]; % Explicit relative path
end

% --- Construct Data Structure ---
DataToSave = struct();
valid_vars = false;

for k = 1:length(varargin)
    % Get the variable name from the caller's workspace
    v_name = inputname(k+1);
    
    if isempty(v_name)
        warning('STAMPS:SaveWarning', 'Argument #%d is an expression, not a named variable. Skipping.', k+1);
        continue;
    end
    
    DataToSave.(v_name) = varargin{k};
    valid_vars = true;
end

if ~valid_vars
    warning('STAMPS:NoVars', 'No valid named variables to save.');
    return;
end

% --- Save Operation ---

if exist(save_name, 'file') == 2
    % === APPEND MODE ===
    try % if raw file is v7.3
        save(save_name, '-struct', 'DataToSave', '-append');
    catch ME
        % if raw file is v7
        if contains(ME.identifier, 'IncompatibleSave') || contains(ME.message, '-v7.3')
            fprintf('Warning: File %s is in legacy format. Upgrading to v7.3 (Large File Support)...\n', save_name);
            
            backup_name = [save_name, '_bak'];
            
            try
                OldData = load(save_name); 
                movefile(save_name, backup_name);
                try
                    % save old data in v7.3
                    save(save_name, '-struct', 'OldData', '-v7.3');
                    % append new data
                    save(save_name, '-struct', 'DataToSave', '-append');
                    
                    delete(backup_name); % clear backup if success
                    fprintf('Success: File upgraded to v7.3. Temporary backup deleted.\n');
                catch SaveError
                    if exist(backup_name, 'file')
                        movefile(backup_name, save_name);
                        fprintf(2, 'Error during save. Restored original file from backup.\n');
                    end
                    rethrow(SaveError); 
                end
                
            catch ConversionError
                error('STAMPS:UpgradeFailed', ...
                      'Failed to load file %s for upgrade. Likely Out of Memory. Error: %s', ...
                      save_name, ConversionError.message);
            end
        else
            rethrow(ME);
        end
    end
else
    % Enforces '-v7.3' (HDF5) for all NEW files.
    save(save_name, '-struct', 'DataToSave', '-v7.3');
end

end