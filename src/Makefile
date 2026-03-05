#################################################################
# StaMPS-HPC Core C/C++ Modules Makefile                        #
#                                                               #
# Commands:                                                     #
#   make [all]      - Compile all core C programs               #
#   make <prog>     - Compile a specific program (e.g., calamp) #
#   make clean      - Remove compiled binaries from src/        #
#   make install    - Compile and copy binaries to INSTALL_DIR  #
#   make uninstall  - Remove binaries from INSTALL_DIR          #
#                                                               #
# Optimized for modern HPC architectures (OpenMP, SIMD)         #
# Author: Mingjia Li                                            #
# Date:   March 2026                                            #
#################################################################

# Compiler settings
CC = g++

# Optimization levels based on specific module requirements
CFLAGS_BASE = -O3 -Wall
CFLAGS_NATIVE = $(CFLAGS_BASE) -march=native
CFLAGS_FASTMATH = $(CFLAGS_NATIVE) -ffast-math

OMPFLAGS = -fopenmp
LDLIBS = -lm

# Target installation directory
INSTALL_DIR = ../bin

# Executable targets
PROGS = calamp pscphase pscdem psclonlat pscheading selpsc_patch selsbc_patch

# Default target
all: $(PROGS)

# 1. calamp: OpenMP + Native optimization
calamp: calamp_opt.c
	$(CC) $(CFLAGS_NATIVE) $(OMPFLAGS) -o $@ $^

# 2. Patch selection: Native optimization + Fast-Math + Math Library (-lm)
selpsc_patch: selpsc_patch_opt.c
	$(CC) $(CFLAGS_FASTMATH) -o $@ $^ $(LDLIBS)

selsbc_patch: selsbc_patch_opt.c
	$(CC) $(CFLAGS_FASTMATH) -o $@ $^ $(LDLIBS)

# 3. pscphase: Base optimization only (matching your specific command)
pscphase: pscphase_opt.c
	$(CC) $(CFLAGS_BASE) -o $@ $^

# 4. Others: Native optimization
pscdem: pscdem_opt.c
	$(CC) $(CFLAGS_NATIVE) -o $@ $^

psclonlat: psclonlat_opt.c
	$(CC) $(CFLAGS_NATIVE) -o $@ $^

pscheading: pscheading.c
	$(CC) $(CFLAGS_NATIVE) -o $@ $^

# Clean compiled binaries from src/
clean:
	@echo "Cleaning up compiled binaries in src..."
	rm -f $(PROGS)

# Install binaries to INSTALL_DIR
install: $(PROGS)
	@if [ ! -d $(INSTALL_DIR) ]; then \
		mkdir -p $(INSTALL_DIR); \
		echo "Created directory: $(INSTALL_DIR)"; \
	fi
	@echo "Installing binaries to $(INSTALL_DIR)..."
	@mv -f $(PROGS) $(INSTALL_DIR)/
	@echo "Installation complete."

# Remove installed binaries from INSTALL_DIR
uninstall:
	@echo "Removing installed binaries from $(INSTALL_DIR)..."
	@cd $(INSTALL_DIR) && rm -f $(PROGS)
	@echo "Uninstallation complete."

.PHONY: all clean install uninstall