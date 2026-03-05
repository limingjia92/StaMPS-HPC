function monitor_func = hpc_log_progress(total_items, step_percent, prefix_str)
%HPC_LOG_PROGRESS Generates a progress monitor for parfor loops.
%   Returns a function handle to be used with parallel.pool.DataQueue.
%
%   Usage:
%       q = parallel.pool.DataQueue;
%       afterEach(q, hpc_log_progress(N, 10, 'Processing'));
%
%   Inputs:
%       total_items  : Total number of iterations
%       step_percent : Output frequency (e.g., 10 for every 10%)
%       prefix_str   : (Optional) Custom prefix for the log message
%
%   Author: Mingjia Li, February 2026
%

    if nargin < 2; step_percent = 10; end
    if nargin < 3; prefix_str = 'Parallel Progress'; end

    % Initialize state variables (Protected inside this closure)
    count = 0;
    last_printed_pct = 0;
    
    % Return the handle to the nested function
    monitor_func = @update_progress;

    % --- Nested Function (Holds the state) ---
    function update_progress(~)
        count = count + 1;
        
        % Calculate current percentage
        pct = floor((count / total_items) * 100);
        
        % Check if we need to print
        if pct >= last_printed_pct + step_percent || count == total_items
            % Update the threshold
            last_printed_pct = pct;
            
            % Print logically (avoiding spam)
            % Using \r could animate it locally, but \n is better for HPC logs
            fprintf('      [%s]: %3d%% complete (%d/%d)\n', ...
                prefix_str, pct, count, total_items);
        end
    end
end