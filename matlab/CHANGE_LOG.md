# StaMPS-HPC MATLAB Processing Optimization Log

**Version:** 1.0.0  
**Status:** High-Performance Hybrid Build  
**Author:** Mingjia Li  
**Focus:** Bottleneck Elimination, Hybrid MEX/OpenMP Acceleration, and I/O Optimization

---

## 1. Executive Summary

This release focuses on resolving the critical computational bottlenecks in the StaMPS InSAR processing chain. By identifying the most time-consuming steps (Step 2 and Step 3) and refactoring them with **OpenMP-accelerated C-MEX kernels** and **Batch Processing strategies**, we have achieved an order-of-magnitude performance improvement.

**Performance Highlight (Tested on Reference Dataset):**
* **Step 3 (PS Selection):** Execution time reduced from **2564s to ~116s** (~22x Speedup).
* **Step 2 (Gamma Estimation):** Execution time reduced from **474s to ~327s** (~1.45x Speedup).

---

## 2. Critical Performance Overhauls (Major Changes)


### [Step 2] Gamma Estimation (`ps_est_gamma_quick.m`)
**Status:** **Core Computation Optimization**
* **The Bottleneck:** Heavy usage of FFT and iterative convolution operations inside the estimation loop.
* **The Solution:** Offloaded heavy math to C-MEX and optimized MATLAB data structures.
* **Key Optimizations:**
    * **MEX Integration:** Replaced `clap_filt.m` with `clap_filt_mex.c` (OpenMP optimized).
    * **Vectorization:** Replaced slow `squeeze` and dimensional permutation operations with direct linear indexing.
    * **I/O Efficiency:** Moved the time-consuming `stamps_save` operation **outside** the main iteration loop, significantly reducing disk I/O latency.

### [Step 3] PS Pixel Selection (`ps_select.m`)
**Status:** **Complete Architectural Refactoring**
* **The Bottleneck:** The original code iterated through candidate pixels sequentially in MATLAB, causing massive interpreter overhead.
* **The Solution:** Implemented a **Batch Processing Strategy**. Instead of looping, the entire pixel grid is passed to a high-performance C-MEX kernel (`clap_filt_patch_mex`) for parallel filtering.
* **Key Optimizations:**
    * **Loop Elimination:** Removed the primary `for` loop over PS candidates.
    * **Hybrid Kernel:** Integrated `clap_filt_patch_mex.c` (OpenMP) for parallel phase filtering.
    * **Optional Approximation:** Added `use_fast_topofit` parameter.
        * `'y'` (Default for HPC): Uses `ps_topofit_mex` (Fastest, ~1.4% numerical divergence).
        * `'n'`: Uses legacy MATLAB topofit (Exact precision, slightly slower).

---

## 3. Stability & Bug Fixes (Minor Changes)

### [Step 1] Initial Data Loading
* **Scripts:** `ps_load_initial_gamma.m`, `sb_load_initial_gamma.m`
* **Changes:** Applied patches to fix stability issues during the initial loading of large-scale binary files.

### [Step 4] Adjacency Weeding
* **Script:** `ps_weed.m`
* **Changes:** * Optimized memory fragmentation during the weeding process.
    * Removed redundant variable allocations to lower peak RAM usage.

### [Step 5] Phase Correction
* **Script:** `ps_correct_phase.m`
* **Changes:** * Fixed critical bugs related to patch merging logic.
    * Streamlined the preparation phase for unwrapping (efficiency patches).

### [Step merge] Patch Merging (`ps_merge_patches.m`)
**Status:** **Architectural Refactoring & Parallelization**
* **The Bottleneck:** The original code processed data "Patch-by-Patch", loading entire workspaces into RAM and appending to global arrays. This caused massive CPU overhead from dynamic resizing and memory thrashing (16,000+ CPU seconds).
* **The Solution:** Implemented a **Variable-Centric Parallel Architecture**.
    1. **Phase 1 (Meta-Scan):** Quickly pre-calculates global indices and sort order without loading heavy data.
    2. **Phase 2 (Parallel Stream):** Processes one variable (e.g., `ph`, `inc`) at a time. Inside each variable task, patches are loaded and processed in parallel using `parfor`.
* **Key Technique:** utilized **Cell Arrays** as intermediate buffers to resolve MATLAB `parfor` variable classification errors and maintain memory stability.
* **Performance:** * **Wall Clock:** **~2.8x Faster** (29m $\to$ 10m).
    * **CPU Efficiency:** **~11.6x Improvement** in effective computation (User Time).
    * **Memory:** Maintained consistent footprint (~15.8GB) matching the original serial version.
    
### [Auxiliary] Interferogram Noise Estimation (`ps_calc_ifg_std.m`)
**Status:** **Mathematical Optimization & Memory Efficiency**
* **The Bottleneck:** The original implementation relied on computationally expensive complex exponential arithmetic (`exp(-j*...)`) and `angle()` functions to compute phase residuals. This generated massive temporary complex arrays, stressing memory bandwidth and CPU FPU.
* **The Solution:** Refactored the core logic to operate entirely in the **Real Number Domain**.
* **Key Optimizations:**
    * **Direct Phase Arithmetic:** Replaced complex multiplications with simple subtraction of phase components (`phase - correction`).
    * **Fast Wrapping:** Implemented a lightweight `mod(x, 2*pi)` logic to replace the slower `angle(complex)` function calls.
    * **Memory Reduction:** Eliminated the creation of massive intermediate complex matrices.
* **Performance:** * **Execution Time:** Reduced from **171s to 41s** (~4.2x Speedup).
    * **Accuracy:** Validated against the original version with negligible deviation (Max Diff < 1e-5 degrees).
---

## 4. Benchmark Results

| Processing Step | Original Time | Optimized Time | Speedup | Note |
| :--- | :--- | :--- | :--- | :--- |
| **Step 2 (Gamma Est)** | 474s | **327s** | **~1.45x** | I/O & Iteration bound |
| **Step 3 (PS Select)** | 2564s | **116s** | **~22.1x** | CPU & Loop bound |

> **Note:** Benchmarks performed on a workstation with SSD storage. Performance gains in Step 3 are primarily due to the elimination of MATLAB loop overheads via MEX batch processing.