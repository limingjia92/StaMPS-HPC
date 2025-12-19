# StaMPS-HPC: High-Performance Parallelized StaMPS

[![Status](https://img.shields.io/badge/Status-Work--in--Progress-orange)](https://github.com/)
[![License](https://img.shields.io/badge/License-GPL--v3.0-blue)](LICENSE)

**StaMPS-HPC** is a performance-optimized derivative of the Stanford Method for Persistent Scatterers (StaMPS). This project aims to refactor the core computational bottlenecks of StaMPS using **OpenMP-accelerated C-MEX** and **Advanced I/O strategies**, specifically designed for processing massive InSAR datasets (e.g., >50 million points).

---

## 🚀 Key Features

* **Parallel Computing:** Multi-threaded execution via OpenMP for both C-binaries and MATLAB-MEX modules.
* **I/O Optimization:** Implemented Block-I/O and Full-Memory reading to minimize disk latency.
* **Hybrid Acceleration:** Replaced heavy MATLAB bottlenecks with high-performance C-MEX kernels.
* **Resource Efficiency:** Optimized memory allocation for large-scale patch merging and phase correction.

## 📊 Performance Benchmark (Preliminary)

| Module | Optimization | Speedup | Note |
| :--- | :--- | :--- | :--- |
| `calamp` | OpenMP + Block I/O | ~70% Faster | Significant on SSD |
| `ps_est_gamma` | C-MEX Refactoring | High | Hybrid Parallelism |
| `ps_select` | Parallel MEX Kernel | ~40% Faster | Configurable Accuracy |

---

## 🛠 Project Structure

* `src/`: Optimized C source code with parallelization. See [src/CHANGE_LOG.md](src/CHANGE_LOG.md).
* `matlab/`: Enhanced MATLAB scripts and MEX interfaces. See [matlab/CHANGE_LOG.md](matlab/CHANGE_LOG.md).
* `bin/`: (Optional) Pre-compiled binaries for supported platforms.

## ⚙️ Quick Start (Development Status)

> **Note:** This project is currently under active development. Some experimental features may require manual configuration.

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/YourUsername/StaMPS-HPC.git](https://github.com/YourUsername/StaMPS-HPC.git)
    ```
2.  **Compile C Binaries:** Follow instructions in `src/CHANGE_LOG.md`.
3.  **Compile MEX Modules:** Follow instructions in `matlab/CHANGE_LOG.md`.

---

## 📑 Intellectual Property & License

### Patent & Copyright
This software is developed as part of a **National Major R&D Project (重点研发项目)**. 
- The parallelization architecture and dynamic I/O scheduling strategies are subject to pending patent applications.
- Software Copyright (软著) is registered for the HPC-optimized implementation.

### Open Source License
This project is licensed under the **GNU General Public License (GPL)**. It is a derivative work based on the original StaMPS software developed by Andy Hooper et al. Users must comply with the original license terms.

---

## ✉️ Contact & Acknowledgement
* **Original Author:** Andy Hooper et al. (Stanford University, University of Leeds).
* **Optimization Work:** Mingjia Li.
* For questions regarding the HPC implementation, please open an **Issue** or contact [Your Email].