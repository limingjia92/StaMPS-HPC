# StaMPS-HPC: Core MATLAB Processing Engine (`matlab`)

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Focus:** Bottleneck Elimination, Hybrid MEX/OpenMP Acceleration, and I/O Optimization 

---

## Module Overview

This directory contains the central MATLAB processing engine for the **StaMPS-HPC** framework. 

Through a rigorous, ground-up refactoring effort, the legacy StaMPS codebase has been aggressively streamlined—**reducing the number of executable `.m` files from over 120 down to 59 essential, highly optimized functions**. 

We have eliminated deprecated MATLAB syntax, eradicated redundant I/O operations, and replaced computationally suffocating loops with heavily vectorized matrix operations and direct C-MEX integrations. The result is a modernized, production-ready InSAR time-series processor capable of handling massive regional datasets.

---

## Performance Benchmarks

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

## Key Engineering Philosophies

1. **Vectorization over Iteration:** O(N) `for`-loops across millions of PS candidates were systematically eradicated and replaced with native MATLAB vectorization (e.g., `accumarray`, `discretize`, and logical masking).
2. **Mathematical Optimization & Memory Efficiency:** Legacy implementations relied on computationally expensive complex exponential arithmetic (`exp(-j*...)`) and `angle()` functions, which generated massive temporary complex arrays and stressed memory bandwidth. These were refactored to operate entirely in the **Real Number Domain** using direct phase subtraction and lightweight `mod(x, 2*pi)` wrappers.
3. **C-MEX & OpenMP Integration:** The most intensive grid searches, spatial convolutions, and filtering loops are automatically handed off to the optimized C-MEX binaries compiled in the `matlab_mex` directory.
4. **Variable-Centric Parallel Architecture:** Process-heavy modules (like patch merging) abandoned the slow "Patch-by-Patch" memory-thrashing approach. Data is now processed one variable at a time using `parfor` and Cell Array Buffering, achieving ~11.6x improvement in CPU efficiency while capping memory footprints.
5. **I/O Consolidation:** Time-consuming `stamps_save` operations were moved outside main iteration loops, drastically reducing disk I/O latency.

---
## Core Processing Chain (Step-by-Step)

### Master Control & Step 1: Initial Ingestion (`stamps`, `load_initial_*`)
* [cite_start]**Workflow & Batch Control (`stamps.m`)**: Decoupled distributed patch processing to natively support parallel job arrays[cite: 10]. [cite_start]Enforced optimized estimation algorithms and integrated external patch list files for precise batch control[cite: 11, 12].
* [cite_start]**Vectorized I/O & Interpolation**: In data loading (GAMMA/ISCE), iterative loops were eradicated in favor of matrix `fread` and fast Regex parsing[cite: 22, 27]. [cite_start]Replaced legacy `griddata` with the highly efficient `interp2` for regular grid interpolation[cite: 23, 28]. [cite_start]Consolidated dependencies (e.g., merging `readparm`) to drastically reduce file handle overhead[cite: 14, 21].

### Step 2: Phase Noise Estimation (`ps_est_gamma_quick`)
* [cite_start]**Hybrid Kernels**: Offloaded heavy FFT and convolution operations to OpenMP-accelerated C-MEX binaries (`clap_filt_mex`)[cite: 31].
* [cite_start]**Matrix Indexing & I/O**: Replaced slow iterative `squeeze` operations with direct linear indexing[cite: 32]. [cite_start]Moved the time-consuming `stamps_save` operation **outside** the main iteration loop, significantly reducing disk I/O latency while maintaining strict numerical consistency[cite: 33, 34].

### Step 3 & 4: Candidate Selection & Weeding (`ps_select`, `ps_weed`)
* [cite_start]**Batch Processing Strategy**: Replaced massive pixel-wise MATLAB loops with a single, highly optimized MEX call (`clap_filt_patch_mex`), eliminating interpreter overhead[cite: 36]. [cite_start]Integrated OpenMP multi-threading for phase filtering[cite: 37].
* [cite_start]**Topofit Toggle**: Introduced `use_fast_topofit` parameter, offering a choice between ultra-fast MEX topofit and exact-precision MATLAB topofit[cite: 37, 38, 39]. 

### Step 5 & 5.5: Phase Correction, Patch Merging & Noise Estimation
* [cite_start]**Two-Phase Parallel Merging (`ps_merge_patches`)**: Implemented a "Meta-Scan -> Parallel Stream" architecture[cite: 41, 42]. [cite_start]Utilized a **Variable-Centric Parfor** approach with Cell Array buffering to resolve dynamic slicing issues and eliminate memory fragmentation[cite: 43, 44, 45]. [cite_start]Selective I/O only loads necessary variables, bypassing entire workspace loads[cite: 46].
* [cite_start]**Real-Domain Arithmetic (`ps_calc_ifg_std`)**: Refactored the noise estimation logic to operate entirely in the Real Number domain (`phase - correction`), dropping expensive complex exponential arithmetic (`exp(-j*...)`)[cite: 48]. Implemented a lightweight `mod(x, 2*pi)` wrapper to replace `angle(complex)`, drastically cutting temporary RAM usage[cite: 49, 50].

### Step 6: 3D Phase Unwrapping (`ps_unwrap` & Sub-modules)
* [cite_start]**Architecture & OS Independence**: Refactored spaghetti code into clean blocks and dynamically parameterized `n_trial_wraps` across sensors[cite: 51, 52]. [cite_start]Eradicated legacy external system calls (e.g., `triangle`) and `uw_nosnaphu`, enforcing a robust Snaphu-based workflow using MATLAB's native `delaunayTriangulation` and `nearestNeighbor`[cite: 53, 60, 61].
* **Parallel Execution & Memory**: Utilized `parfor` for concurrent interferogram gridding and statistical cost generation[cite: 57, 70]. Leveraged `snaphu -S` (Tile Mode) with `NPROC=1` for single-threaded efficiency within workers[cite: 71]. Eliminated massive history matrices (~99% RAM saving) during smooth unwrapping via JIT-inlined objective functions[cite: 65, 66].
* [cite_start]**Advanced Solver & Vectorization**: Replaced pixel-wise loops with vectorized matrix operations (`accumarray`)[cite: 55, 58]. [cite_start]Replaced the iterative `lscov` solver with a pre-calculated Generalized Least Squares (GLS) operator $H$ and implemented block chunking to minimize RAM footprints[cite: 73, 74].

### Step 7: Error Estimation (SCLA, Orbit & Atmosphere)
* **L1/L2 Vectorization (`ps_calc_scla`)**: 
    * [cite_start]**L2 Mode:** Replaced iterative `lscov` with a vectorized GLS matrix operator[cite: 77]. 
    * [cite_start]**L1 Mode:** Completely eradicated the computationally suffocating `fminsearch`, replacing it with a fully vectorized Iteratively Reweighted Least Squares (IRLS) algorithm to boost speed and prevent local minima traps[cite: 78].
* **Deramping (`ps_deramp`)**: Implemented a "Common Valid Pixel" strategy for consistent orbital plane estimation[cite: 75]. [cite_start]Replaced iterative solvers with vectorized QR decomposition (`\`) using BLAS Level 3[cite: 75].
* [cite_start]**Graph Traversal & Projection (`ps_smooth_scla`)**: Replaced loop-heavy neighborhood searches with `accumarray` (dropping time to seconds)[cite: 81]. [cite_start]Substituted giant matrix transpose divisions with a pre-calculated projection operator ($Data * H_{op}'$) to eliminate OOM risks[cite: 82]. 

### Step 8: Spatial Filtering & Kriging (`ps_scn_filt`, `ps_scn_filt_krig`)
* [cite_start]**Topology & KD-Tree Search**: Shifted to in-memory `delaunayTriangulation` and robust KD-Tree spatial indexing (`rangesearch`, `knnsearch`), boosting speed and fixing heuristic neighborhood bugs[cite: 85, 87, 89, 92].
* **Custom Micro-Solvers**: Bypassed `lsqcurvefit` entirely. Built a custom, inlined Levenberg-Marquardt solver with analytical Jacobians for variogram fitting (10x speedup)[cite: 91, 99]. Replaced ill-conditioned normal equations with QR decomposition (`A\y`) to fix catastrophic $2\pi$ unwrapping errors[cite: 90].
* [cite_start]**Vectorized Kriging (`ps_kriging`)**: Suppressed OOM crashes via implicit expansion (dropping `repmat`), enforced $N_{max}$ data truncation, and hardcoded theoretical variogram models directly into loops to strip evaluation overheads[cite: 93, 94, 95].

---

## Post-Processing, Models & Visualization

### Atmospheric & Ionospheric Corrections (`sb_invert_tca`, `sb_invert_iono`)
* [cite_start]**String-Driven Architecture**: Refactored legacy `if-elseif` chains into performant `switch-case` logic using readable string flags (e.g., `'a_era5'`, `'i_split'`), while maintaining integer backward compatibility[cite: 106, 108, 109, 116].
* **OOM Prevention**: Replaced double-transpose `lscov` operations with a memory-efficient projection operator ($H_{op}$)[cite: 105, 112]. Extended native support to ERA5, Mcandis stacking, and Split-Spectrum/TEC map integrations[cite: 110, 115].

### High-Performance Interactive Rendering (`ps_plot`, `ts_plotdiff`, `profile_plot`)
* **Zero I/O & Memory Encapsulation**: Legacy workspace pollution was eliminated. [cite_start]Time-Series UIs now utilize figure-bound `appdata` for encapsulated, conflict-free interactions[cite: 141, 158]. [cite_start]Functions like `ps_plot_ifg` and `env_oscilator_corr` pass pre-loaded structures directly, averting secondary disk loads and memory bloat[cite: 143, 154].
* **Vectorized UI Engines**: Eradicated outdated $O(N)$ pixel dilation and swath binning loops, replacing them with ultra-fast `accumarray` and `discretize` functions for instantaneous visualization[cite: 144, 166].
* [cite_start]**Advanced UX**: Implemented robust global marker tracking to seamlessly manage crosshairs across tools[cite: 157, 162]. [cite_start]Modernized rendering palettes via `cptcmap` for robust GMT `.cpt` parsing and integrated interactive dual-axis swath profiles[cite: 142, 146, 152].---
