# StaMPS-HPC Source Code Optimization Log

**Version:** 1.0.0  
**Status:** Performance-Enhanced Build  
**Author:** Mingjia Li  
**Focus:** I/O Throughput, Memory Management, and OpenMP Parallelization

---

## 1. Amplitude Calibration (`calamp`)

**Optimization Strategy:**
* **Parallelization:** Refactored the core processing loop using **OpenMP** to enable multi-threaded execution.
* **I/O Architecture:** Transitioned from legacy line-by-line reading to a **Full-Memory / Block I/O** strategy. This minimizes the frequency of disk system calls by loading large chunks of `.rslc` data directly into RAM.
* **Overhead Reduction:** Significantly reduced the CPU wait time associated with disk latency.

**Performance Impact:**
* **Speedup:** Achieved approximately **70% reduction** in execution time.
* **Hardware Recommendation:** As the process is I/O-bound, utilizing NVMe SSDs instead of traditional HDDs is highly recommended for maximum throughput.

**Compilation & Usage:**
Bash
# Compilation (Requires OpenMP)
g++ -o calamp calamp_opt.c -O3 -fopenmp -march=native

# Execution (Example with 4 threads)
export OMP_NUM_THREADS=4 
./calamp calamp.in 67546 calamp.out s 1 > calamp.log

## 2. PS/SB Selection (selsbc_patch & selpsc_patch)

**Optimization Strategy:**
* **Buffered Reading:** Implemented a LineInBufferMax mechanism to handle data in optimized memory blocks rather than continuous small streams.
* **Memory Efficiency:** Optimized the allocation of internal structures to better handle large-scale datasets.

**Performance Impact:**
* **Speedup:** Achieved approximately 40% reduction in execution time.
* **Execution Note:** Due to the high I/O saturation, running multiple concurrent instances on a single HDD may lead to disk contention and sub-optimal performance.

**Compilation & Usage:**
Bash

g++ -O3 -march=native -ffast-math -o selsbc_patch selsbc_patch_opt.c -lm
g++ -O3 -march=native -ffast-math -o selpsc_patch selpsc_patch_opt.c -lm

## 3. Phase Stability Estimation (pscphase)

**Optimization Strategy:**
* **Seek vs. Read Trade-off:** Modified the I/O logic to prioritize sequential reading over random seeking.
* **Latent Seek Reduction:** By increasing the file system input volume (FS Inputs), the program reduces mechanical disk head movement (Seek latency).
* **Tunable Parameter:** Introduced READ_VS_SEEK_THRESHOLD (Default: 2MB). The system automatically decides whether to "read through" or "seek" based on the gap between data points.

**Performance Impact:**
* **Benefit:** Most effective for high-density PS point sets. While the raw speedup is marginal on modern SSDs, it significantly improves efficiency on enterprise-level HDDs.

**Compilation & Usage:**

Bash

g++ -O3 -o pscphase pscphase_opt.c

## 4. DEM and Coordinates Search (pscdem & psclonlat)
**Optimization Strategy:**
* **Algorithm Refactoring:** Updated search algorithms to improve lookup efficiency during the coordinate assignment phase.
* **Instruction Optimization:** Applied -march=native to leverage modern CPU instruction sets.

**Performance Impact:**
* **Result:* Provides a cleaner execution profile, though the absolute time reduction is minor as these modules are not the primary bottleneck in the StaMPS workflow.

**Compilation & Usage:**

Bash

g++ -O3 -march=native -o pscdem pscdem_opt.c
g++ -O3 -march=native -o psclonlat psclonlat_opt.c