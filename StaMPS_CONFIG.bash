#!/bin/bash
#################################################################
# StaMPS-HPC Environment Configuration Script                   #
#                                                               #
# Purpose: Adds StaMPS-HPC directories to PATH and MATLABPATH.  #
# Usage:   source StaMPS_CONFIG.bash                            #
#                                                               #
# Author:  Mingjia Li                                           #
# Date:    March 2026                                           #
#################################################################

# Dynamically get the absolute path of the StaMPS-HPC root directory
STAMPS_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Add the bin directory to the system PATH for C executables
export PATH="$STAMPS_ROOT/bin:$PATH"

# Add matlab and matlab_mex directories to MATLABPATH
export MATLABPATH="$STAMPS_ROOT/matlab:$STAMPS_ROOT/matlab_mex:$MATLABPATH"