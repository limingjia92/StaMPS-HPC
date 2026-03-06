# StaMPS-HPC: High-Performance Parallelized InSAR Time-Series Framework

[![Status](https://img.shields.io/badge/Status-Active-brightgreen)](https://github.com/limingjia92/StaMPS-HPC)
[![License](https://img.shields.io/badge/License-GPL--v3.0-blue)](LICENSE)

**StaMPS-HPC** is a performance-optimized derivative of the Stanford Method for Persistent Scatterers (StaMPS). This project aims to refactor the core computational bottlenecks of the original StaMPS using OpenMP-accelerated C-MEX and Advanced I/O strategies, specifically designed for processing massive InSAR datasets (e.g., >50 million points). 

Through a rigorous, ground-up refactoring effort across C/C++, MATLAB, and Bash architectures, we have transformed the legacy codebase into a modernized, production-ready time-series processor capable of handling regional-scale InSAR analysis with unprecedented efficiency.

---

## Key Features & Innovations

* **Hybrid Parallel Computing**: Integrated **OpenMP multi-threading** across C/C++ binaries and MATLAB C-MEX modules, offloading computationally heavy tasks (e.g., FFTs, spatial convolutions) directly to compiled C kernels.
* **Extreme I/O & Memory Control**: Replaced legacy line-by-line file reading with Block-I/O strategies to minimize disk latency. Solved MATLAB memory thrashing during massive patch merging via a novel "Variable-Centric" `parfor` architecture.
* **Algorithmic Refactoring**: Streamlined the MATLAB engine from 120+ to 59 core scripts. Eradicated slow O(N) loops in favor of native vectorization (`accumarray`, `discretize`) and refactored suffocating complex exponential arithmetic into the Real Number Domain.
* **Modernized Ecosystem & UX**: Enabled multi-core concurrent patch extraction for ISCE2 and GAMMA. Overhauled the visualization engine with 30+ standard GMT Colormaps (`.cpt`) and introduced new interactive utilities (`ps_export_kmz`, `profile_plot`).
---

## Performance Benchmarks

Significant engineering focus was placed on eliminating algorithmic bottlenecks present throughout the entire processing chain. Benchmarks demonstrate massive performance gains achieved through architectural refactoring.

### C/C++ Core Binaries
*(Tested on SSD environments)*
* **Amplitude Calibration (`calamp`)**: ~70% reduction in execution time due to Block I/O and OpenMP integration.
* **PS/SB Selection (`selpsc` & `selsbc`)**: ~40% reduction in execution time via buffered reading logic.

### MATLAB Processing Engine
*(Tested on Reference Dataset with HDD storage)*
| Processing Step | Original Time | HPC Time | Speedup | Core Optimization |
| :--- | :--- | :--- | :--- | :--- |
| **Step 3 (`ps_select`)** | 2564s | **116s** | **22.1x** | Loop Elimination, C-MEX Hybrid Kernel |
| **Step 5.5 (`ps_calc_ifg_std`)** | 171s | **41s** | **4.2x** | Real-Domain Math, Fast Wrapping |
| **Step 5.5 (`ps_merge_patches`)**| 1470s | **600s** | **2.8x** | Variable-Centric Parfor, Cell Buffering |
| **Step 6 (`ps_unwrap`)** | 1068s | **411s** | **2.6x** | Parallel Execution, Snaphu Optimization |
| **Step 8 (`ps_scn_filt_krig`)** | 3149s | **408s** | **7.7x** | Variogram Fitting, Spatial Kriging KD-Tree |

---

## Project Structure

The project is modularized into five highly optimized sub-components:

* `bin/`: High-performance Bash wrapper scripts serving as the bridge between ISCE2/GAMMA and StaMPS, featuring concurrent patch processing capabilities (`jobmaxnum`).
* `src/`: Optimized C/C++ core source codes responsible for initial candidate selection and geometry extraction.
* `matlab/`: The modernized central MATLAB processing engine, featuring deep algorithmic vectorization, refined mathematical solvers, and an upgraded interactive plotting suite.
* `matlab_mex/`: High-performance C-MEX modules designed to bypass MATLAB interpreter bottlenecks.
* `cptfiles/`: A comprehensive library of GMT Color Palette Tables seamlessly parsed by the internal engine to enhance visualization aesthetics.

---

## Quick Start & Installation

> **Note:** Processing on NVMe SSDs is strongly recommended to maximize throughput and prevent disk contention during multi-core processing.

### 1. Installation & Compilation
Clone the repository and compile all C/C++ core binaries and MATLAB C-MEX modules. The root `Makefile` will automatically trigger the respective build processes in both the `src` and `matlab_mex` directories. (Ensure `gcc` with OpenMP support and `matlab` are available in your system `$PATH`).

``` bash
git clone https://github.com/limingjia92/StaMPS-HPC.git
cd StaMPS-HPC
make install
```

### 2. Environment Configuration
To use the StaMPS-HPC commands from any directory, you need to configure the environment variables.

**Temporary Configuration (Current Terminal Session Only):**
``` bash
source StaMPS_CONFIG.bash
```

**Permanent Configuration:**
If you want to apply these variables permanently across all new terminal sessions, you must add the absolute path of the configuration script to your `~/.bashrc` (or `~/.zshrc`) file. Open your shell profile and append the following line:
``` bash
source /absolute/path/to/StaMPS-HPC/StaMPS_CONFIG.bash
```
*(Do not forget to replace `/absolute/path/to/` with your actual system path, and run `source ~/.bashrc` afterward).*

---

## Intellectual Property & License

### Copyright
Software Copyright is officially registered for this HPC-optimized implementation. 

### Open Source License
This project is licensed under the **GNU General Public License (GPL)**. It is a derivative work based on the original StaMPS software. Users must comply with the original license terms.

---

## Acknowledgement & Contact

* **Optimization & Refactoring:** Mingjia Li
* **Original Authors:** Andy Hooper et al. (Stanford University, University of Leeds)
* **Acknowledgements:** Special thanks to **Jianbao Sun** and **Kejie Chen** for their invaluable support, guidance, and contributions to the underlying concepts of this optimization project.

For questions regarding the HPC implementation, parallelization configuration, or bug reports, please open an **Issue** on GitHub.