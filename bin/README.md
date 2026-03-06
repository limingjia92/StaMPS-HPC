# StaMPS-HPC: Pre-processing Wrapper Scripts (`bin`)

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Date:** March 2026

---

## Module Overview

This directory contains the essential Bash shell scripts that serve as the high-performance bridge between modern InSAR processors (e.g., **ISCE2**, **GAMMA**) and the StaMPS-HPC core processing chain.

These wrapper scripts have been extensively refactored to support **parallel patch processing**, optimize data I/O, and seamlessly invoke the optimized C/C++ core binaries (e.g., `calamp`, `pscdem`, `pscphase`) compiled in the `src` directory.

---

## Detailed Script Log

### 1. ISCE2 Stack Preparation (`make_isce_stack_ps.sh` & `make_isce_stack_sbas.sh`)
* **Function:** Prepares the Single Reference (PS) or Small Baselines (SBAS) stacks from ISCE2 SLC outputs (typically processed via `stackSentinel.py`).
* **Key Optimization:** * Completely replaces the legacy `make_single_reference_stack_isce` and `make_small_baselines_isce` scripts.
    * Streamlines the generation of VRT/XML headers and logical symlinks without duplicating massive SLC datasets.
    * Incorporates automatic interferogram generation (`imageMath.py`) and quicklook visualization options.

### 2. StaMPS Data Extraction (`prep_stamps_isce.sh` & `prep_stamps_gamma.sh`)
* **Function:** Extracts amplitude dispersion, phase, height, and geographic coordinates for PS/SB candidates, organizing the data into `PATCH_*` directories for subsequent MATLAB processing.
* **Key Optimization:**
    * **Parallelization:** Introduced the `jobmaxnum` parameter, enabling concurrent patching of massive datasets across multiple CPU cores.
    * **Direct Binary Integration:** Automatically routes data to the newly optimized C/C++ binaries, stripping away legacy dependencies.
    * **GAMMA Specifics:** Added optional heading and incidence angle extraction (`pscheading`) utilizing look vector results.

---

## Workflow & Usage

The general workflow for processing data using these scripts involves stacking the SLCs first, followed by patch preparation. 

### Example Workflow for ISCE2 (PS Processing):
```bash
# Step 1: Prepare the ISCE2 stack (Reference date: 20220104)
make_isce_stack_ps.sh 20220104 merged/SLC merged/geom_reference merged/baselines

# Step 2: Extract candidates and prepare StaMPS patches with parallelization 
# (e.g., Da_thresh: 0.4, Patches: 4x8, Parallel Jobs: 5)
prep_stamps_isce.sh INSAR_20220104 work_dir 0.4 4 8 5
```

(For detailed parameter inputs, simply run any of the scripts without arguments to view the usage menu).

---

## Performance & Configuration

### I/O Bottleneck Warning

While the prep_stamps_*.sh scripts fully support multi-core parallelization via the jobmaxnum argument, performance during extraction is heavily bottlenecked by disk read/write speeds.

* **NVMe SSDs:** Highly recommended. You can safely increase jobmaxnum to match your CPU cores for massive speedups.

* **Mechanical HDDs:** Increasing jobmaxnum may decrease overall performance due to severe disk thrashing. Limit parallel jobs to 1 or 2 if operating on HDDs.

### Path Configuration

Ensure that this bin directory is accessible in your system's $PATH. If you have executed the top-level StaMPS_CONFIG.bash script, this requirement is already handled automatically.