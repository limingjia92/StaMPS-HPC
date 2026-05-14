# Changelog

All notable changes to the StaMPS-HPC project will be documented in this file.

## [Unreleased] - 2026-03-25

### Added
- **`matlab/ps_weed.m`**: Introduced an Auto-Fallback mechanism for `psver` management (2026-03-24).  
  - **Purpose**: Prevents "file not found" crashes when re-running the script with modified parameters 
- **`matlab/plot_baselines.m`**: Introduced a unified and modernized spatiotemporal baseline network plotting tool (2026-03-25).  
  - **Purpose**: Combines legacy PS and SBAS plotting scripts, and automated high-resolution PNG exports
- **`bin/make_isce_stack_ps.sh` & `bin/make_isce_stack_sbas.sh`**: Added automated generation of high-resolution PNG previews for Amplitude and Phase interferograms (2026-03-25).
  - **Purpose**: Enables rapid visual quality inspection immediately after stack generation. 
- **`matlab/ps_plot.m`**: Major data export overhaul and UX/UI enhancements (2026-03-26).
  - **Purpose**: Enhanced user experience during data extraction by:
    1) Adding 'c'/'csb' flags for direct cumulative deformation extraction in millimeters (mm).
    2) Dynamically generating physically variable names (e.g., `disp_velocity`) and saving `units` for background exports (`plot_flag = -1`).
    3) Resolving the ambiguous `mean_v.mat` output by explicitly saving `velocity` and `lscov_coeffs` instead of raw mathematical variables.
- **`matlab/ps_load_initial_isce.m` & `matlab/sb_load_initial_isce.m`**: Introduced automatic calculation and export of Look Angle (`la1.mat`) derived from the incidence angle and satellite orbit altitude (2026-04-02).
  - **Purpose**: Bridges the data gap between ISCE2 outputs and downstream atmospheric correction tools (like TRAIN) that strictly require Look Angle geometry rather than Incidence Angle.
- **`bin/make_isce_stack_ps.sh`, `bin/make_isce_stack_sbas.sh`, `bin/prep_stamps_isce.sh`, `matlab/ps_parms_initial.m`, & `matlab/sb_parms_initial.m`**: Implemented an end-to-end pipeline for the automated extraction, propagation, and saving of the `UTC_sat` variables (2026-04-02).
  - **Purpose**: Eliminates manual metadata lookups by automatically passing the satellite overpass time from the raw ISCE files directly into the StaMPS `parms.mat` structure, seamlessly preparing the dataset for time-dependent atmospheric correction models (e.g., ERA5 or GACOS).
- **`matlab/ps_gacos2tca.m`**: Introduced a standalone utility to natively assimilate GACOS tropospheric delay products (GeoTIFF/Binary) into the processing workflow (2026-04-02).
  - **Purpose**: Computes the TCA phase matrix for both PS and SBAS modes via direct geographic interpolation and strict reference point alignment, entirely eliminating the reliance on heavy external toolboxes like TRAIN.
  
### Changed
- **`matlab/uw_stat_costs.m`**: Removed the experimental Pre-Unwrapping Filtering feature (2026-03-25).  
  - **Purpose**: Ensures 100% mathematical consistency with original StaMPS cost calculations.
- **`matlab/ps_merge_patches.m`**: Completely refactored using a "Streamed Pre-allocation & Global Mapping" architecture (2026-04-02).
  - **Purpose**: Eliminates `parfor` and `vertcat` operations to guarantee 100% OOM prevention on massive datasets (e.g., >30M points), locking peak memory to the theoretical minimum.
- **`matlab/ps_calc_ifg_std.m`**: Redesigned with a "Row-Block Streaming" (Dimensional Inversion) architecture (2026-04-02).
  - **Purpose**: Resolves severe HDF5 "I/O Thrashing" by sequentially reading data row-by-row. Strictly caps peak RAM under 10GB while maximizing multi-threaded CPU vectorization via implicit expansion.
- **`matlab/ps_unwrap.m`, `matlab/uw_3d.m` & `matlab/uw_grid_wrapped.m`**: Completely overhauled the phase unwrapping I/O pipeline using "Row-Block Streaming" and a "Two-Phase Gridding" architecture (2026-04-03).
  - **Purpose**: Eradicates OOM crashes and disk I/O deadlocks on massive datasets (>30M points) by bypassing full-matrix memory loading and optimizing `matfile` data access for `parfor`.
- **`matlab/uw_sb_unwrap_space_time.m` & `matlab/uw_stat_costs.m`**: Upgraded temporal smoothing and statistical cost generation with memory-safe `parfor` slicing and implicit expansion (2026-04-03).
  - **Purpose**: Accelerates 3D unwrapping and SNAPHU preparation while strictly preventing parfor variable broadcast spikes and `repmat` memory overhead.
- **`matlab/sb_invert_uw.m`**: Redesigned SBAS inversion with "Two-Pass Block-IO" covariance calculation and streamed matrix multiplication (2026-04-03).
  - **Purpose**: Eliminates >130GB RAM requirements during full-matrix loading and accelerates computation using a vectorized Generalized Least Squares (GLS) operator.
- **`matlab/ps_calc_scla.m` & `matlab/ps_smooth_scla.m`**: Overhauled the SCLA estimation and smoothing pipelines with "Streamed Casting", Block-IO Projection, and implicit expansion (2026-04-04).
  - **Purpose**: Eradicates >140GB RAM spikes by strictly confining memory-heavy operations (e.g., `double(diff())`, `repmat`, and pseudo-inverse projection) into chunked, low-memory footprints. Massively accelerates solvers via vectorized GLS/IRLS and native Delaunay graph traversal (`accumarray`).
- **`matlab/ps_deramp.m`**: Implemented "Strict Precision Casting" and vectorized QR decomposition for orbital ramp estimation (2026-04-04).
  - **Purpose**: Prevents fatal OOM crashes (e.g., 65GB+ spikes on 30M points) by restricting `double` precision strictly to the mathematical solver and immediately downgrading reconstructed ramps to `single`, preventing implicit type expansion during phase subtraction.
- **`matlab/ps_scn_filt.m` & `matlab/ps_scn_filt_krig.m`**: Overhauled the SCN filtering pipeline utilizing "Node-Based Convolution" and "Parfor Loop Inversion" (2026-04-04).
  - **Purpose**: Eradicates catastrophic OOM crashes (>380GB RAM spikes & >768GB parfor broadcasts) by filtering lightweight nodes instead of massive edge matrices, and enforcing strict memory slicing in spatial Kriging.
- **`matlab/ps_plot.m`**: Eradicated >150GB RAM spikes during velocity inversion and bootstrapping by implementing Block-IO GLS solvers and strict precision bounds (2026-04-04).
  - **Purpose**: Prevents fatal OOM crashes 
- **`matlab/ts_plot.m` & `matlab/ts_plotdiff.m`**: Overhauled time-series interactive extraction tools with "Lazy Load Memory Mapping" and "Appdata Binding" (2026-04-06).
  - **Purpose**: Eradicates UI freezing and OOM crashes when extracting data from massive (e.g., >29GB) TS files by utilizing `matfile` pointers with contiguous indexing bypasses. Replaces fragile global variables with figure-bound appdata, implements strict SBAS-network dimensional safeguards, and visually integrates the zero-deformation master date.
- **`matlab/profile_plot.m`**: Completely refactored the interactive swath profile tool with "Vectorized Binning" and "Dual-Axis Integration" (2026-04-06).
  - **Purpose**: Eliminates legacy `O(N)` looping bottlenecks by using `accumarray` and `discretize` for instantaneous swath aggregation. Drastically reduces memory footprints by fetching lightweight 1D velocity arrays directly from disk, and modernizes the UX with dynamic dual-Y-axis elevation rendering and global crosshair tracking.
  
### Fixed
- **`matlab_mex/Makefile`**: Added auto-detection for the `-batch` execution flag (2026-03-23).  
  - **Purpose**: Restores seamless MEX compilation compatibility for legacy MATLAB versions (R2018b and earlier) . 
- **`matlab/ps_merge_patches.m`**: Replaced hardcoded `n_cols_ifg` with dynamic column size acquisition for complex variables (2026-03-24).  
  - **Purpose**: Resolves "Matrix dimensions must agree" errors during patch resampling for standard PS datasets. 
- **`matlab/ps_scn_filt_krig.m`**: Removed the `krig_atmo` conditional logic and enforced the Kriging APS filtering workflow unconditionally (2026-05-13)
  - **Purpose**: Keeps the optimized implementation consistent with the original ps_scn_filt_krig behavior, where krig_atmo was hardcoded to true, and avoids accidental execution of the untested non-Kriging fallback branch.
- **`matlab/ps_parms_default.m`**:: Removed the default `krig_atmo` parameter setting. (2026-05-13)
  - **Purpose**: Eliminates an obsolete configuration option that no longer affects ps_scn_filt_krig, reducing ambiguity and preventing users from assuming that non-Kriging APS filtering is supported in this workflow.
- **`matlab/ps_scn_filt_krig.m`**: Added robust spatial Kriging safeguards and final output validation (2026-05-14).
  - **Purpose**: Prevents all-NaN SCN columns caused by unrealistically small variogram ranges or sparse Kriging neighborhoods by enforcing geometry-based range lower bounds, nearest-neighbor fallback, robust SCN re-referencing, and pre-save hard checks for invalid `ph_scn_slave` / `ph_noise_res` outputs.
