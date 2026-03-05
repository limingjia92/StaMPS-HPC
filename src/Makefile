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
#################################################################

# Compiler settings
CC = g++
CFLAGS = -O3 -march=native -Wall
OMPFLAGS = -fopenmp

# Target installation directory
INSTALL_DIR = ../bin

# Executable targets
PROGS = calamp pscphase pscdem psclonlat pscheading selpsc_patch selsbc_patch

# Default target
all: $(PROGS)

# Compile with OpenMP support
calamp: calamp.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $^

# Compile with SIMD/native optimization
pscphase: pscphase_opt.c
	$(CC) $(CFLAGS) -o $@ $^

pscdem: pscdem_opt.c
	$(CC) $(CFLAGS) -o $@ $^

psclonlat: psclonlat_opt.c
	$(CC) $(CFLAGS) -o $@ $^

pscheading: pscheading.c
	$(CC) $(CFLAGS) -o $@ $^

# Compile patch selection modules
selpsc_patch: selpsc_patch.c
	$(CC) $(CFLAGS) -o $@ $^

selsbc_patch: selsbc_patch.c
	$(CC) $(CFLAGS) -o $@ $^

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
	@cp -f $(PROGS) $(INSTALL_DIR)/
	@echo "Installation complete."

# Remove installed binaries from INSTALL_DIR
uninstall:
	@echo "Removing installed binaries from $(INSTALL_DIR)..."
	@cd $(INSTALL_DIR) && rm -f $(PROGS)
	@echo "Uninstallation complete."

.PHONY: all clean install uninstall