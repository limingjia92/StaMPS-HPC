function ps_parms_initial()
%PS_PARMS_INITIAL Initialize parms for PS processing (HPC Standardized)
%   Sets small_baseline_flag='n' and detects processor type.
%   Andy Hooper, Jan 2008
%   Updated January 2026

    parmfile = 'parms.mat';
    
    % 1. Initialize Structure
    parms = struct('Created', date);
    parms.small_baseline_flag = 'n'; % PS mode (Single Reference)

    % 2. Processor Detection 
    processor_file = 'processor.txt';
    proc_path = '';
    
    if exist(processor_file, 'file')
        proc_path = processor_file;
    elseif exist(['..' filesep processor_file], 'file')
        proc_path = ['..' filesep processor_file];
    end
    
    if ~isempty(proc_path)
        raw_text = fileread(proc_path);
        parms.insar_processor = strtrim(raw_text);
    else
        % Default for HPC environment
        parms.insar_processor = 'isce'; 
        fprintf('Info: processor.txt not found. Defaulting to ISCE.\n');
    end

    % 3. Save to Disk
    try
        save(parmfile, '-struct', 'parms');
    catch ME
        error('Stamps:WriteError', 'Could not save %s. Check permissions.\n%s', parmfile, ME.message);
    end

    % 4. Apply Defaults
    ps_parms_default;

end