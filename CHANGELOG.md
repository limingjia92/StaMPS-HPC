# Changelog

All notable changes to the StaMPS-HPC project will be documented in this file.

## [Unreleased] - 2026-03-25

### Added
- **`matlab/ps_weed.m`**: Introduced an Auto-Fallback mechanism for `psver` management. Prevents "file not found" crashes when re-running the script with modified parameters (2026-03-24). 

### Changed
- **`matlab/uw_stat_costs.m`**: Removed the experimental Pre-Unwrapping Filtering feature. Ensures 100% mathematical consistency with original StaMPS cost calculations (2026-03-25). 

### Fixed
- **`matlab/ps_merge_patches.m`**: Replaced hardcoded `n_cols_ifg` with dynamic column size acquisition for complex variables. Resolves "Matrix dimensions must agree" errors during patch resampling for standard PS datasets. (2026-03-24). 
- **`matlab_mex/Makefile`**: Added auto-detection for the `-batch` execution flag.  Restores seamless MEX compilation compatibility for legacy MATLAB versions (R2018b and earlier) (2026-03-23). 