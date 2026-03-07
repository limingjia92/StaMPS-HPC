function build_mex()
% BUILD_MEX Automatically compile C-MEX modules for StaMPS-HPC.
%
% Purpose: Compiles the optimized C source files into MATLAB executable (MEX)
%          binaries with OpenMP acceleration and high-level optimizations.
%          This script is designed to be called in a headless MATLAB session.
%
% Author:  Mingjia Li
% Date:    March 2026
% -------------------------------------------------------------------------

% Suppress MATLAB-level warnings (specifically the GCC version mismatch warning)
warning('off', 'all');

disp('===================================================');
disp('   Starting StaMPS-HPC MEX Compilation...          ');
disp('===================================================');

% Define compilation flags (Note: double single-quotes are used to escape inside eval)
flags = '-R2018a CFLAGS=''$CFLAGS -fopenmp -O3'' LDFLAGS=''$LDFLAGS -fopenmp''';

try
    disp('-> Compiling clap_filt_mex.c ...');
    eval(['mex ', flags, ' clap_filt_mex.c']);
    
    disp('-> Compiling ps_topofit_mex.c ...');
    eval(['mex ', flags, ' ps_topofit_mex.c']);
    
    disp('-> Compiling clap_filt_patch_mex.c ...');
    eval(['mex ', flags, ' clap_filt_patch_mex.c']);

    disp('-> Compiling smooth_arcs_mex.c ...');
    eval(['mex ', flags, ' smooth_arcs_mex.c']);
    
    disp('===================================================');
    disp('   Compilation Successful!                         ');
    disp('===================================================');
catch ME
    disp('===================================================');
    disp('   Compilation Failed!                             ');
    disp(ME.message);
    disp('===================================================');
    % Throw an error code to the OS so the Makefile knows the build failed
    exit(1); 
end

end