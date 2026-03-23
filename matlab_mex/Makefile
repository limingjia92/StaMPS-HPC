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

# Auto-detect MATLAB version: Check if '-batch' is in the help documentation
# Count matching lines using grep -c. > 0 means R2019a or newer; 0 means older versions
SUPPORT_BATCH := $(shell matlab -help 2>/dev/null | grep -c "\-batch")

# Dynamically set commands and arguments based on the version
ifeq ($(SUPPORT_BATCH), 0)
	# For R2018b and older: Use legacy mode, requires an explicit 'exit' to prevent process hanging
	MATLAB_CMD = matlab -nojvm -nodisplay -nosplash -nodesktop -r
	RUN_ARGS = "build_mex; exit"
else
	# For R2019a and newer: Use native non-interactive clean batch mode
	MATLAB_CMD = matlab -nojvm -nodisplay -batch
	RUN_ARGS = "build_mex"
endif

# Default target
all:
	@echo "Launching headless MATLAB to compile MEX files..."
	@$(MATLAB_CMD) $(RUN_ARGS)

# Clean target
clean:
	@echo "Cleaning up compiled MEX binaries..."
	rm -f *.mexa64 *.mexw64 *.o

.PHONY: all clean