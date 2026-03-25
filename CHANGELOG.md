# Changelog

All notable changes to the StaMPS-HPC project will be documented in this file.

## [Unreleased] - 2026-03-25

### Added
- **`matlab/ps_weed.m`**: Introduced an Auto-Fallback mechanism for `psver` management (2026-03-24).  
    Prevents "file not found" crashes when re-running the script with modified parameters 

### Changed
- **`matlab/uw_stat_costs.m`**: Removed the experimental Pre-Unwrapping Filtering feature (2026-03-25). 
    Ensures 100% mathematical consistency with original StaMPS cost calculations.

### Fixed
- **`matlab_mex/Makefile`**: Added auto-detection for the `-batch` execution flag (2026-03-23).  
    Restores seamless MEX compilation compatibility for legacy MATLAB versions (R2018b and earlier) . 
- **`matlab/ps_merge_patches.m`**: Replaced hardcoded `n_cols_ifg` with dynamic column size acquisition for complex variables (2026-03-24). 
    Resolves "Matrix dimensions must agree" errors during patch resampling for standard PS datasets. 
    
