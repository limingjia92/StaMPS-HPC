#################################################################
# StaMPS-HPC MATLAB MEX Modules Makefile (Headless Mode)        #
#                                                               #
# Commands:                                                     #
#   make [all]      - Launch headless MATLAB to compile MEX     #
#   make clean      - Remove compiled MEX binaries              #
#                                                               #
# Author: Mingjia Li                                            #
# Date:   March 2026                                            #
#################################################################

# Use '-batch' for clean, non-interactive execution (R2019a+)
MATLAB_CMD = matlab -nojvm -nodisplay -batch

# Default target
all:
	@echo "Launching headless MATLAB to compile MEX files..."
	@$(MATLAB_CMD) "build_mex"

# Clean target
clean:
	@echo "Cleaning up compiled MEX binaries..."
	rm -f *.mexa64 *.mexw64 *.o

.PHONY: all clean