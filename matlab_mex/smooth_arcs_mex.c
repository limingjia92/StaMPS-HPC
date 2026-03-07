/*
 * ==============================================================================
 * File: smooth_arcs_mex.c
 * Original Author: Andy Hooper (ps_weed.m)
 * Optimization Author: Mingjia Li
 * Date: March 2026
 * License: GNU General Public License (GPL)
 * ==============================================================================
 * DESCRIPTION:
 * This is a high-performance C-MEX implementation of the arc phase smoothing and 
 * noise estimation logic from 'ps_weed.m'. It processes millions of arcs across 
 * temporal baselines to filter out noise prior to standard deviation thresholding.
 *
 * OPTIMIZATION STRATEGY:
 * 1. Loop Fusion & Cache Optimization: Eliminated massive temporary matrices 
 * (Memory Bandwidth Bound) by fusing phase multiplication, angle extraction, 
 * and least-squares fitting into a single cache-friendly loop.
 * 2. Hardcoded WLS Solver: Replaced MATLAB's computationally expensive 'lscov' 
 * function with a hardcoded Cramer's Rule solver for 2x2 weighted linear 
 * systems, drastically reducing interpreter overhead.
 * 3. Phase Wrapping Optimization: Replaced expensive 'angle(exp(1j * x))' calls 
 * with a pure arithmetic modulo-based phase wrapping algorithm (fmod).
 * 4. OpenMP Parallelization: Parallelized the processing of millions of arcs 
 * across multi-core CPUs using shared memory, entirely avoiding the OOM 
 * risks associated with MATLAB's 'parfor'.
 * * COMPILATION:
 * mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" smooth_arcs_mex.c
 * ==============================================================================
 */

#include "mex.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* 1. Input Validation */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("StaMPS:smooth_arcs_mex:nrhs", "Requires exactly 4 inputs: dph_space, dph_mean_all, day_ifg, W.");
    }

    /* 2. Extract pointers and dimensions from MATLAB inputs */
    mxComplexDouble *dph_space = mxGetComplexDoubles(prhs[0]);
    size_t n_edge = mxGetM(prhs[0]);
    size_t n_use = mxGetN(prhs[0]);

    mxComplexDouble *dph_mean_all = mxGetComplexDoubles(prhs[1]);
    double *day_ifg = mxGetDoubles(prhs[2]);
    double *W = mxGetDoubles(prhs[3]);

    /* 3. Initialize Output Matrix (Single Precision Complex) */
    plhs[0] = mxCreateUninitNumericMatrix(n_edge, n_use, mxSINGLE_CLASS, mxCOMPLEX);
    mxComplexSingle *dph_smooth = mxGetComplexSingles(plhs[0]);

    /* 4. Determine OpenMP Thread Count for Workspace Allocation */
    int max_threads = 1;
    #ifdef _OPENMP
    max_threads = omp_get_max_threads();
    #endif
    
    /* 5. Allocate memory manually (C99 requirement) */
    /* Thread-local workspace array to prevent race conditions during OpenMP execution */
    double *y_workspace = (double*)malloc(max_threads * n_use * sizeof(double));
    double *t = (double*)malloc(n_use * sizeof(double));
    double *w = (double*)malloc(n_use * sizeof(double));

    /* Check for successful memory allocation */
    if (y_workspace == NULL || t == NULL || w == NULL) {
        if (y_workspace) free(y_workspace);
        if (t) free(t);
        if (w) free(w);
        mexErrMsgIdAndTxt("StaMPS:smooth_arcs_mex:memory", "Failed to allocate memory for workspaces.");
    }

    /* 6. Outer Loop: Iterate through each interferogram (i1) sequentially */
    for (size_t i1 = 0; i1 < n_use; ++i1) {
        double Sw = 0, Swt = 0, Swtt = 0;
        
        /* Pre-compute temporal baselines (t) and WLS normal equation components */
        for (size_t k = 0; k < n_use; ++k) {
            t[k] = day_ifg[i1] - day_ifg[k];
            w[k] = W[k + i1 * n_use]; /* W is column-major in MATLAB */
            Sw += w[k];
            Swt += w[k] * t[k];
            Swtt += w[k] * t[k] * t[k];
        }
        
        /* Determinant for 2x2 Weighted Least Squares (Cramer's Rule) */
        double det = Sw * Swtt - Swt * Swt;
        int valid_det = fabs(det) > 1e-15 ? 1 : 0;

        /* 7. Inner Loop: Process all Delaunay arcs in parallel */
        #pragma omp parallel for
        for (size_t e = 0; e < n_edge; ++e) {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            double* y = &y_workspace[tid * n_use]; 

            /* Construct complex mean phase for the current edge and interferogram */
            double complex mean_all_e_i1 = dph_mean_all[e + i1 * n_edge].real + 
                                           dph_mean_all[e + i1 * n_edge].imag * I;

            /* =========================================================
               Phase 1: First Weighted Least Squares Regression
               ========================================================= */
            double Swy = 0, Swty = 0;
            for (size_t k = 0; k < n_use; ++k) {
                size_t idx = e + k * n_edge;
                double complex space_e_k = dph_space[idx].real + dph_space[idx].imag * I;
                
                /* Multiply by conjugate of the weighted mean */
                double complex prod = space_e_k * conj(mean_all_e_i1);
                
                /* Extract phase angle (equivalent to MATLAB's angle) */
                y[k] = carg(prod); 
                
                Swy += w[k] * y[k];
                Swty += w[k] * t[k] * y[k];
            }

            /* Solve 2x2 system analytically */
            double m0 = 0, m1 = 0;
            if (valid_det) {
                m0 = (Swtt * Swy - Swt * Swty) / det;
                m1 = (Sw * Swty - Swt * Swy) / det;
            }

            /* =========================================================
               Phase 2: Phase Wrapping and Second Regression
               ========================================================= */
            double Swy_new = 0, Swty_new = 0;
            for (size_t k = 0; k < n_use; ++k) {
                /* Calculate residual phase from the first linear fit */
                double y_new = y[k] - (m0 + m1 * t[k]);
                
                /* Fast arithmetic phase wrapping to [-PI, PI] */
                y_new = fmod(y_new + M_PI, 2.0 * M_PI);
                if (y_new < 0) y_new += 2.0 * M_PI;
                y_new -= M_PI;
                
                Swy_new += w[k] * y_new;
                Swty_new += w[k] * t[k] * y_new;
            }

            /* Solve second 2x2 system (we only need the intercept m0_new) */
            double m0_new = 0;
            if (valid_det) {
                m0_new = (Swtt * Swy_new - Swt * Swty_new) / det;
            }

            /* =========================================================
               Phase 3: Reconstruction and Output
               ========================================================= */
            double total_phase = m0 + m0_new;
            
            /* Reconstruct the smoothed complex value manually */
            double complex out_val = mean_all_e_i1 * (cos(total_phase) + sin(total_phase) * I);

            /* Write back to the pre-allocated MATLAB output matrix */
            size_t out_idx = e + i1 * n_edge;
            dph_smooth[out_idx].real = (float)creal(out_val);
            dph_smooth[out_idx].imag = (float)cimag(out_val);
        }
    }

    /* 8. Free allocated memory to prevent memory leaks */
    free(y_workspace);
    free(t);
    free(w);
}