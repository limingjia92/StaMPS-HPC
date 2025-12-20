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

---

## 4. Benchmark Results

| Processing Step | Original Time | Optimized Time | Speedup | Note |
| :--- | :--- | :--- | :--- | :--- |
| **Step 2 (Gamma Est)** | 474s | **327s** | **~1.45x** | I/O & Iteration bound |
| **Step 3 (PS Select)** | 2564s | **116s** | **~22.1x** | CPU & Loop bound |

> **Note:** Benchmarks performed on a workstation with SSD storage. Performance gains in Step 3 are primarily due to the elimination of MATLAB loop overheads via MEX batch processing.