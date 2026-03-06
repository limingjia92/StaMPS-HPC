# StaMPS-HPC: Core MATLAB Processing Engine (`matlab`)

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Focus:** Bottleneck Elimination, Hybrid MEX/OpenMP Acceleration, and I/O Optimization

---

## Module Overview

This directory contains the central MATLAB processing engine for the **StaMPS-HPC** framework. 

Through a rigorous, ground-up refactoring effort, the legacy StaMPS codebase has been aggressively streamlined—reducing the number of executable `.m` files from over 120 down to 59 essential, highly optimized functions. 

We have eliminated deprecated MATLAB syntax, eradicated redundant I/O operations, and replaced computationally suffocating loops with heavily vectorized matrix operations and direct C-MEX integrations. The result is a modernized, production-ready InSAR time-series processor capable of handling massive regional datasets.

---

## Performance Benchmarks

Significant engineering focus was placed on eliminating algorithmic bottlenecks present throughout the entire processing chain. Benchmarks demonstrate massive performance gains achieved through architectural refactoring, OpenMP-accelerated C-MEX kernels, and Batch Processing strategies.

*(Tested on Testing Dataset with HDD storage)*

| Processing Module | Original Time | HPC Time | Speedup | Core Optimization Strategy |
| :--- | :--- | :--- | :--- | :--- |
| **Step 2** (`ps_est_gamma_quick`) | 474s | **327s** | **1.45x** | MEX Integration, Vectorization, I/O Efficiency |
| **Step 3** (`ps_select`) | 2564s | **116s** | **22.1x** | Loop Elimination, Hybrid Kernel |
| **Step 5.5** (`ps_calc_ifg_std`) | 171s | **41s** | **4.2x** | Direct Phase Arithmetic, Fast Wrapping, Memory Reduction |
| **Step 5.5** (`ps_merge_patches`) | 1470s | **600s** | **2.8x** | Variable-Centric Parfor, Cell Array Buffering |
| **Step 6** (`ps_unwrap`) | 1068s | **411s** | **2.6x** | Structural Optimization, Parallel Execution, Snaphu Optimization |
| **Step 7** (`ps_calc_scla`) | 29s | **14s** | **2.1x** | L2 Vectorization & Precision, L1 Vectorized IRLS |
| **Step 7** (`ps_smooth_scla`) | 252s | **47s** | **5.4x** | Topology Generation, Vectorized Graph Traversal, Memory Projection |
| **Step 8** (`ps_scn_filt`) | 529s | **80s** | **6.6x** | Topology Generation, Vectorized Filtering & Solving |
| **Step 8** (`ps_scn_filt_krig`) | 3149s | **408s** | **7.7x** | Numerical Stability, Variogram Fitting, Spatial Kriging |

---

## Key Engineering Philosophies

1. **Vectorization over Iteration:** O(N) `for`-loops across millions of PS candidates were systematically eradicated and replaced with native MATLAB vectorization (e.g., `accumarray`, `discretize`, and logical masking).
2. **Mathematical Optimization & Memory Efficiency:** Legacy implementations relied on computationally expensive complex exponential arithmetic (`exp(-j*...)`) and `angle()` functions, which generated massive temporary complex arrays. These were refactored to operate entirely in the Real Number Domain using direct phase subtraction and lightweight `mod(x, 2*pi)` wrappers.
3. **C-MEX & OpenMP Integration:** Intensive grid searches, spatial convolutions, and filtering loops are automatically handed off to optimized C-MEX binaries compiled in the `matlab_mex` directory.
4. **Variable-Centric Parallel Architecture:** Process-heavy modules abandoned the slow "Patch-by-Patch" memory-thrashing approach. Data is now processed one variable at a time using `parfor` and Cell Array Buffering, achieving massive CPU efficiency while capping memory footprints.
5. **I/O Consolidation:** Time-consuming `stamps_save` operations were moved outside main iteration loops, drastically reducing disk I/O latency.

---

## Core Processing Chain (Step-by-Step)

### Master Control & Step 1: Initial Ingestion (`stamps`, `load_initial_*`)
* **Workflow & Batch Control (`stamps.m`)**: Decoupled distributed patch processing to natively support parallel job arrays. Enforced optimized estimation algorithms and integrated external patch list files for precise batch control.
* **Vectorized I/O & Interpolation**: In data loading (GAMMA/ISCE), iterative loops were eradicated in favor of matrix `fread` and fast Regex parsing. Replaced legacy `griddata` with the highly efficient `interp2` for regular grid interpolation. Consolidated dependencies (e.g., merging `readparm`) to drastically reduce file handle overhead.

### Step 2: Phase Noise Estimation (`ps_est_gamma_quick`)
* **Hybrid Kernels**: Offloaded heavy FFT and convolution operations to OpenMP-accelerated C-MEX binaries (`clap_filt_mex`).
* **Matrix Indexing & I/O**: Replaced slow iterative `squeeze` operations with direct linear indexing. Moved the time-consuming `stamps_save` operation outside the main iteration loop, significantly reducing disk I/O latency while maintaining strict numerical consistency.

### Step 3 & 4: Candidate Selection & Weeding (`ps_select`, `ps_weed`)
* **Batch Processing Strategy**: Replaced massive pixel-wise MATLAB loops with a single, highly optimized MEX call (`clap_filt_patch_mex`), eliminating interpreter overhead. Integrated OpenMP multi-threading for phase filtering.
* **Topofit Toggle**: Introduced `use_fast_topofit` parameter, offering a choice between ultra-fast MEX topofit and exact-precision MATLAB topofit.

### Step 5 & 5.5: Phase Correction, Patch Merging & Noise Estimation
* **Two-Phase Parallel Merging (`ps_merge_patches`)**: Implemented a "Meta-Scan -> Parallel Stream" architecture. Utilized a Variable-Centric Parfor approach with Cell Array buffering to resolve dynamic slicing issues and eliminate memory fragmentation. Selective I/O only loads necessary variables.
* **Real-Domain Arithmetic (`ps_calc_ifg_std`)**: Refactored the noise estimation logic to operate entirely in the Real Number domain (`phase - correction`), dropping expensive complex exponential arithmetic (`exp(-j*...)`). Implemented a lightweight `mod(x, 2*pi)` wrapper to replace `angle(complex)`, drastically cutting temporary RAM usage.

### Step 6: 3D Phase Unwrapping (`ps_unwrap` & Sub-modules)
* **Architecture & OS Independence**: Refactored spaghetti code into clean blocks and dynamically parameterized `n_trial_wraps` across sensors. Eradicated legacy external system calls (e.g., `triangle`) and `uw_nosnaphu`, enforcing a robust Snaphu-based workflow using MATLAB's native `delaunayTriangulation` and `nearestNeighbor`.
* **Parallel Execution & Memory**: Utilized `parfor` for concurrent interferogram gridding and statistical cost generation. Leveraged `snaphu -S` (Tile Mode) with `NPROC=1` for single-threaded efficiency within workers. Eliminated massive history matrices (~99% RAM saving) during smooth unwrapping via JIT-inlined objective functions.
* **Advanced Solver & Vectorization**: Replaced pixel-wise loops with vectorized matrix operations (`accumarray`). Replaced the iterative `lscov` solver with a pre-calculated Generalized Least Squares (GLS) operator ($H$) and implemented block chunking to minimize RAM footprints.

### Step 7: Error Estimation (SCLA, Orbit & Atmosphere)
* **L1/L2 Vectorization (`ps_calc_scla`)**: 
    * **L2 Mode:** Replaced iterative `lscov` with a vectorized GLS matrix operator. 
    * **L1 Mode:** Completely eradicated the computationally suffocating `fminsearch`, replacing it with a fully vectorized Iteratively Reweighted Least Squares (IRLS) algorithm to boost speed and prevent local minima traps.
* **Deramping (`ps_deramp`)**: Implemented a "Common Valid Pixel" strategy for consistent orbital plane estimation. Replaced iterative solvers with vectorized QR decomposition (`\`) using BLAS Level 3.
* **Graph Traversal & Projection (`ps_smooth_scla`)**: Replaced loop-heavy neighborhood searches with `accumarray` (dropping time to seconds). Substituted giant matrix transpose divisions with a pre-calculated projection operator ($Data * H_{op}'$) to eliminate OOM risks. 

### Step 8: Spatial Filtering & Kriging (`ps_scn_filt`, `ps_scn_filt_krig`)
* **Topology & KD-Tree Search**: Shifted to in-memory `delaunayTriangulation` and robust KD-Tree spatial indexing (`rangesearch`, `knnsearch`), boosting speed and fixing heuristic neighborhood bugs.
* **Custom Micro-Solvers**: Bypassed `lsqcurvefit` entirely. Built a custom, inlined Levenberg-Marquardt solver with analytical Jacobians for variogram fitting (10x speedup). Replaced ill-conditioned normal equations with QR decomposition (`A\y`) to fix catastrophic unwrapping errors.
* **Vectorized Kriging (`ps_kriging`)**: Suppressed OOM crashes via implicit expansion (dropping `repmat`), enforced $N_{max}$ data truncation, and hardcoded theoretical variogram models directly into loops to strip evaluation overheads.

---

## Post-Processing, Models & Visualization

### Atmospheric & Ionospheric Corrections (`sb_invert_tca`, `sb_invert_iono`)
* **String-Driven Architecture**: Refactored legacy `if-elseif` chains into performant `switch-case` logic using readable string flags (e.g., `'a_era5'`, `'i_split'`), while maintaining integer backward compatibility.
* **OOM Prevention**: Replaced double-transpose `lscov` operations with a memory-efficient projection operator ($H_{op}$). Extended native support to ERA5, Mcandis stacking, and Split-Spectrum/TEC map integrations.

### High-Performance Interactive Rendering (`ps_plot`, `ts_plotdiff`, `profile_plot`)
* **Zero I/O & Memory Encapsulation**: Legacy workspace pollution was eliminated. Time-Series UIs now utilize figure-bound `appdata` for encapsulated, conflict-free interactions. Functions like `ps_plot_ifg` and `env_oscilator_corr` pass pre-loaded structures directly, averting secondary disk loads and memory bloat.
* **Vectorized UI Engines**: Eradicated outdated $O(N)$ pixel dilation and swath binning loops, replacing them with ultra-fast `accumarray` and `discretize` functions for instantaneous visualization.
* **Advanced UX**: Implemented robust global marker tracking to seamlessly manage crosshairs across tools. Modernized rendering palettes via `cptcmap` for robust GMT `.cpt` parsing and integrated interactive dual-axis swath profiles.

---

## New Additions & Original Enhancements

To further enrich the StaMPS-HPC ecosystem, several entirely new tools and utilities were developed from scratch to improve workflow transparency, large-scale data export, and visual analysis:

* **`hpc_log_progress`**: A custom utility that generates a robust progress monitor specifically designed for `parfor` loops, ensuring visibility and progress tracking during heavy parallel execution.
* **`ps_export_kmz`**: Enables the direct export of large-scale PS point clouds into rasterized KMZ files. This overcomes the severe rendering limitations of standard KMLs when dealing with millions of scatterers on platforms like Google Earth.
* **`plot_patch_kml_gamma` & `plot_patch_kml_isce`**: New utilities to instantly generate KML files that visualize the spatial boundaries and distribution of StaMPS patches derived directly from GAMMA or ISCE workflows.
* **`profile_plot`**: A fully interactive swath profile plotting tool. It features mouse-driven dual-point selection, dynamic line rendering, and synchronized dual-Y axes for topography and deformation velocity. It utilizes `discretize` and `accumarray` for instantaneous $O(1)$ binning, completely eradicating legacy looping delays.