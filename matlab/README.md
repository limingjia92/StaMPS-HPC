# StaMPS-HPC: Core MATLAB Processing Engine (`matlab`)

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Focus:** Bottleneck Elimination, Hybrid MEX/OpenMP Acceleration, and I/O Optimization 

---

## 📌 Module Overview

This directory contains the central MATLAB processing engine for the **StaMPS-HPC** framework. 

Through a rigorous, ground-up refactoring effort, the legacy StaMPS codebase has been aggressively streamlined—**reducing the number of executable `.m` files from over 120 down to 59 essential, highly optimized functions**. 

We have eliminated deprecated MATLAB syntax, eradicated redundant I/O operations, and replaced computationally suffocating loops with heavily vectorized matrix operations and direct C-MEX integrations. The result is a modernized, production-ready InSAR time-series processor capable of handling massive regional datasets.

---

## 🚀 Performance Benchmarks

Significant engineering focus was placed on eliminating the algorithmic bottlenecks present throughout the entire processing chain. The following benchmarks demonstrate the massive performance gains achieved through architectural refactoring, OpenMP-accelerated C-MEX kernels, and Batch Processing strategies.

*(Tested on Reference Dataset with SSD storage)*

| Processing Module | Original Time | HPC Time | Speedup | Core Optimization Strategy |
| :--- | :--- | :--- | :--- | :--- |
| **Step 2** (`ps_est_gamma_quick`) | 474s | **327s** | **1.45x** | MEX Integration, Vectorization, I/O Efficiency |
| **Step 3** (`ps_select`) | 2564s | **116s** | **22.1x** | Loop Elimination, Hybrid Kernel |
| **Step 5.5** (`ps_calc_ifg_std`) | 171s | **41s** | **4.2x** | Direct Phase Arithmetic, Fast Wrapping, Memory Reduction |
| **Step 5.5** (`ps_merge_patches`) | 29m | **10m** | **2.8x** | Variable-Centric Parfor, Cell Array Buffering |
| **Step 6** (`ps_unwrap`) | *** | *** | **\*\*x** | Structural Optimization, Parallel Execution, Snaphu Optimization |
| **Step 7** (`ps_calc_scla`) | 29s | **14s** | **2.1x** | L2 Vectorization & Precision, L1 Vectorized IRLS |
| **Step 7** (`ps_smooth_scla`) | 252s | **47s** | **5.4x** | Topology Generation, Vectorized Graph Traversal, Memory Projection |
| **Step 8** (`ps_scn_filt`) | 529s | **80s** | **6.6x** | Topology Generation, Vectorized Filtering & Solving |
| **Step 8** (`ps_scn_filt_krig`) | 3149s | **408s** | **7.7x** | Numerical Stability, Variogram Fitting, Spatial Kriging |

---

## 🧠 Key Engineering Philosophies

1. **Vectorization over Iteration:** O(N) `for`-loops across millions of PS candidates were systematically eradicated and replaced with native MATLAB vectorization (e.g., `accumarray`, `discretize`, and logical masking).
2. **Mathematical Optimization & Memory Efficiency:** Legacy implementations relied on computationally expensive complex exponential arithmetic (`exp(-j*...)`) and `angle()` functions, which generated massive temporary complex arrays and stressed memory bandwidth. These were refactored to operate entirely in the **Real Number Domain** using direct phase subtraction and lightweight `mod(x, 2*pi)` wrappers.
3. **C-MEX & OpenMP Integration:** The most intensive grid searches, spatial convolutions, and filtering loops are automatically handed off to the optimized C-MEX binaries compiled in the `matlab_mex` directory.
4. **Variable-Centric Parallel Architecture:** Process-heavy modules (like patch merging) abandoned the slow "Patch-by-Patch" memory-thrashing approach. Data is now processed one variable at a time using `parfor` and Cell Array Buffering, achieving ~11.6x improvement in CPU efficiency while capping memory footprints.
5. **I/O Consolidation:** Time-consuming `stamps_save` operations were moved outside main iteration loops, drastically reducing disk I/O latency.

---

## ⚙️ Core Processing Chain (Step-by-Step)

### 🎛️ Master Control & Step 1: Initial Ingestion
* **`stamps.m`**: Workflow restructuring, argument simplification, and batch control. Completely restructured to decouple distributed patch processing, supporting parallel job arrays.
* **`ps_load_initial_*` & `sb_load_initial_*`**: Dependency consolidation and algorithmic efficiency. Unified support for GAMMA and ISCE processors, replacing iterative loops with block `fread` operations and fast Regex parsing.

### 📉 Step 2: Phase Noise Estimation (`ps_est_gamma_quick`)
* **Optimizations:** MEX Integration, Vectorization, I/O Efficiency. Substituted native `gamma` algorithms with scalable sparse-matrix logic. Replaced `clap_filt.m` with `clap_filt_mex.c` (OpenMP optimized).

### 🎯 Step 3: Candidate Selection (`ps_select`)
* **Optimizations:** Loop Elimination, Hybrid Kernel. The entire pixel grid is now passed to a high-performance C-MEX kernel (`clap_filt_patch_mex`) for parallel filtering. 
* **Feature:** Added `use_fast_topofit` parameter to toggle between `ps_topofit_mex` (fastest) and legacy MATLAB topofit (exact precision).

### 🌿 Step 4 & 5: Weeding & Phase Correction
* **`ps_weed` / `ps_correct_phase`:** Small optimizations, memory fragmentation reduction, and bug fixes related to patch merging logic. Replaced notoriously slow built-in `intersect` and `setdiff` commands with logical masking array operations.

### 🔗 Step 5.5: Patch Merging & Noise Estimation
* **`ps_merge_patches`:** Implemented a Two-Phase Architecture (Meta-Scan + Parallel Stream). Utilized Variable-Centric Parfor and Cell Array Buffering for extreme CPU efficiency and memory reduction.
* **`ps_calc_ifg_std`:** Direct Phase Arithmetic, Fast Wrapping, and Memory Reduction.

### 🗺️ Step 6: 3D Phase Unwrapping (`ps_unwrap`)
* **Optimizations:** Structural Optimization, Dynamic Parameterization, Parallel Execution, and Memory Efficiency. Robust 3D_FULL implementation with Pre-Unwrapping Filtering and system-level Snaphu Optimization. Integrated fast spatial search KD-Trees.

### 🌦️ Step 7: Error Estimation (`ps_calc_scla`, `ps_smooth_scla`)
* **Optimizations:** Architecture Refactoring, L2 Vectorization & Precision, L1 Vectorized IRLS. Topology Generation, Vectorized Graph Traversal, and Memory Projection. 

### 🌪️ Step 8: Spatial Filtering (`ps_scn_filt`, `ps_scn_filt_krig`)
* **Optimizations:** Topology Generation, Vectorized Filtering & Solving, Parallel Processing. For kriging, improvements focused on Numerical Stability, Variogram Fitting, and Spatial Kriging.

---

## 📊 Post-Processing, Models & Visualization

* **Framework Control (`ps_parms_default`)**: Framework Optimization, Parameter Cleanup, Annotation.
* **Corrections (`sb_invert_tca`, `sb_invert_iono`)**: Architecture Alignment, String-Driven Logic, Extended Model Support.
* **Interactive Plotting (`ps_plot`, `ts_plotdiff`, `profile_plot`)**: Parameter Parsing, Logic Purification, TS Data Binding, UI Overhaul & Profiling. Features interactive marker systems, GMT standard `.cpt` colormaps, and instantaneous swath binning via `discretize` and `accumarray`.