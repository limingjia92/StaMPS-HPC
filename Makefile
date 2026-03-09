#################################################################
# StaMPS-HPC Top-Level Makefile                                 #
#                                                               #
# Commands:                                                     #
#   make [all]      - Compile C core and MATLAB MEX modules     #
#   make install    - Compile all and deploy C binaries to bin/ #
#   make clean      - Remove all compiled binaries              #
#   make uninstall  - Remove installed binaries from bin/       #
#                                                               #
# Author: Mingjia Li                                            #
# Date:   March 2026                                            #
#################################################################

.PHONY: all install clean uninstall c_core mex_core

# Default target: Compile all modules without deployment
all: c_core mex_core
	@echo ">>> All StaMPS-HPC modules compiled successfully!"

# Install target: Compile MEX and deploy C binaries
install: mex_core
	@echo ">>> Building and installing C Core Modules..."
	@$(MAKE) -C src install
    @echo ">>> Setting executable permissions for Bash scripts..."
	@chmod 755 bin/*.sh
	@echo ">>> StaMPS-HPC installation complete!"

# Compile C source code only
c_core:
	@echo ">>> Compiling C Core Modules..."
	@$(MAKE) -C src all

# Compile MATLAB MEX source code only
mex_core:
	@echo ">>> Compiling MATLAB MEX Modules..."
	@$(MAKE) -C matlab_mex all

# Clean compiled binaries in all subdirectories
clean:
	@echo ">>> Cleaning C Core Modules..."
	@$(MAKE) -C src clean
	@echo ">>> Cleaning MATLAB MEX Modules..."
	@$(MAKE) -C matlab_mex clean
	@echo ">>> All temporary binaries cleaned!"

# Uninstall deployed binaries from the bin directory
uninstall:
	@echo ">>> Uninstalling C Core Modules..."
	@$(MAKE) -C src uninstall
	@echo ">>> Uninstallation complete!"