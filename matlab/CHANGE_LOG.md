# StaMPS-HPC MATLAB Processing Optimization Log

**Version:** 1.0.0  
**Status:** Hybrid MATLAB/C-MEX Performance Build  
**Author:** Mingjia Li  
**Focus:** Algorithmic Refactoring, OpenMP-Accelerated MEX, and Hybrid Workflow Efficiency

---

## 1. Summary of Enhancements

This update represents a major performance overhaul of the StaMPS processing chain (Steps 1-5). By transitioning computationally intensive MATLAB bottlenecks to **OpenMP-accelerated C-MEX** implementations, we have achieved a significant reduction in total processing time while maintaining high numerical stability.

---

## 2. Detailed Step-by-Step Changes

### [Step 1] Initial Data Loading
* **Scripts:** `ps_load_initial_gamma.m`, `sb_load_initial_gamma.m`
* **Changes:** Resolved stability bugs during the Gamma data ingestion phase to ensure robust handling of large-scale datasets.

### [Step 2] Gamma Estimation (Core Computation)
* **Script:** `ps_est_gamma_quick.m`
* **Optimization:** Complete refactoring using hybrid MEX-C integration.
* **Key Replacements:**
    * Replaced MATLAB `clap_filt.m` with high-performance `clap_filt_mex.c`.
    * Replaced MATLAB `ps_topofit.m` with high-performance `ps_topofit_mex.c`.
* **Technology:** Leveraged **OpenMP multi-threading** within the C-MEX kernel to handle pixel-wise calculations in parallel.

### [Step 3] PS Pixel Selection
* **Script:** `ps_select.m`
* **Optimization:** Re-engineered the selection logic to support external C-MEX calls.
* **Key Replacements:**
    * Replaced MATLAB `clap_filt_patch.m` with `clap_filt_patch_mex.c`.
    * Integrated optional `ps_topofit_mex.c` for accelerated topographic fitting.
* **Accuracy Trade-off:** * The C-based `topofit` introduces a minor numerical divergence of **~1.4e-02 (1.4%)** compared to the original MATLAB implementation. 
    * **Toggle:** Users can enable this via `setparm('use_fast_topofit', 'y')`. The default is set to `'n'` for maximum precision.

### [Step 4] Adjacency Weeding
* **Script:** `ps_weed.m`
* **Changes:** Minor algorithmic optimizations to reduce redundant loops and memory fragmentation.

### [Step 5] Phase Correction & Patch Merging
* **Script:** `ps_correct_phase.m`
* **Changes:** Fixed critical bugs and applied efficiency patches to the phase unwrapping preparation and patch merging logic.

---

## 3. Compilation Instructions (MEX)

To activate the acceleration modules, compile the source files within the MATLAB environment.

### Prerequisites
* **Compiler:** GCC 6.3.x or later (Linux/WSL2).
* **Architecture:** MATLAB R2018a or later (64-bit).

### Build Commands
Execute the following in the MATLAB Command Window:

matlab
% Compile Step 2 Modules (clap_filt & ps_topofit)
mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" clap_filt_mex.c
mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" ps_topofit_mex.c

% Compile Step 3 Modules (clap_filt_patch)
mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" clap_filt_patch_mex.c

##4. Usage & Configuration
Once compiled, the resulting .mexa64 (Linux) or .mexw64 (Windows) files must be present in the MATLAB search path. The processing scripts will automatically prioritize the MEX implementations over the original .m files.

Note: Performance gains are most noticeable on workstations with high core counts (e.g., 16+ cores) and datasets exceeding 10 million PS points.