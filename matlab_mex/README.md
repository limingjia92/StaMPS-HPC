# MATLAB MEX Optimization Module

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Date:** March 2026

---

## 📌 Module Overview

This directory contains the **C-MEX (MATLAB Executable)** source codes designed to replace the most computationally intensive MATLAB scripts in the StaMPS processing chain. 

By offloading Fast Fourier Transforms (FFT), convolution, and iterative grid searches to compiled C code, we achieve significant performance gains while maintaining numerical consistency with the original algorithms.

## 🚀 Quick Start: Automated Compilation

We have integrated an automated, headless MATLAB compilation workflow using a unified `Makefile` and the `build_mex.m` function. **You no longer need to open the MATLAB GUI or manually input `mex` commands.**

Ensure the `matlab` command is available in your system's PATH, then simply run the following commands in this directory:

```bash
# Launch headless MATLAB to compile all C-MEX files with OpenMP acceleration
make

# Clean up compiled MEX binaries and temporary object files
make clean
```

Note: The Makefile automatically invokes build_mex.m in batch mode (-batch) and applies the optimal -fopenmp and -O3 flags.

## 📝 Detailed Change Log

### 1. `ps_topofit_mex.c` (Replaces `ps_topofit.m`)
* **Bottleneck Solved:** The original MATLAB script used a slow iterative loop for "Coarse Search" estimation of look-angle errors.
* **Key Optimization:** * Implemented **OpenMP** to parallelize the search over thousands of PS candidates.
    * Forced **Double Precision Accumulators** to fix numerical instability (cycle slips) observed when porting logic from MATLAB to C.
    * Strict adherence to MATLAB's `sum(psdph==0)==0` rejection logic.

### 2. `clap_filt_mex.c` (Replaces `clap_filt.m`)
* **Bottleneck Solved:** Heavy usage of `fft2`, `ifft2`, and `filter2` inside nested loops for sliding window processing.
* **Key Optimization:**.
    * **Trigonometric Lookup Tables (LUT):** Pre-calculated sin/cos values for DFT, eliminating redundant math library calls.
    * **Separable Convolution:** Decomposed the 7x7 Gaussian kernel into 1D filters, reducing operation count by ~70%.
    * Parallelized window processing using OpenMP.

### 3. `clap_filt_patch_mex.c` (Replaces `clap_filt_patch.m`)
* **Bottleneck Solved:** Repeated MEX invocation overhead when processing millions of PS candidates individually.
* **Key Optimization:**
    * **Batch Processing:** The loop over candidates is moved inside the C kernel.
    * **Kernel Normalization:** Added specific normalization logic to ensure the convolution output magnitude matches MATLAB's built-in filter exactly.
    * Added explicit **NaN handling** to safeguard against edge-case data corruption.

## 🚀 Usage & Configuration
### File Priority
Once compiled, the resulting binary files (e.g., .mexa64 on Linux or .mexw64 on Windows) must be present in the MATLAB search path.

The StaMPS processing scripts are designed to automatically prioritize these MEX implementations over the original .m files if they verify the presence of the compiled binaries.

### Performance Note
* **Hardware Recommendation:** Performance gains are most noticeable on workstations with high core counts (e.g., 16+ cores) and when processing datasets exceeding 10 million PS points.