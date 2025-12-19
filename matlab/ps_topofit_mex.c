/*
 * ps_topofit_mex.c
 * Final Stability Version
 *
 * Revisions:
 * 1. Matches MATLAB's 'sum(psdph==0)==0' strict logic (UNAMBIGUOUS).
 * 2. FORCES DOUBLE PRECISION ACCUMULATION in Coarse Search for numerical stability, 
 * even when input is single. This is a compromise to prevent accumulated K_ps
 * errors (cycle slips) which diverge in the iterative MATLAB loop.
 * 3. Maintains Double precision for Trigonometry.
 *
 * Compile: mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" ps_topofit_mex.c
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Coarse Search: DOUBLE Accumulators, Double Trig (Improved Stability) */
void coarse_search_double_accum(
    int n_trials, int n_valid, 
    double *trial_mult, double *bperp_local, double factor,
    double *psdph_r, double *psdph_i, double sum_abs,
    double *max_coh_out, int *best_k_idx_out) 
{
    // max_coh uses double for stability during search
    double max_coh = -1.0; 
    int best_k_idx = 0;
    
    // No need for Temp Float Buffers

    double factor_d = factor;

    for (int t = 0; t < n_trials; t++) {
        double tm = trial_mult[t];
        
        // Accumulators are DOUBLE for maximum stability
        double sr = 0.0; 
        double si = 0.0;
        
        for (int k = 0; k < n_valid; k++) {
            double ang = -(bperp_local[k] * factor_d) * tm;
            
            // Trig and multiplication remain in Double
            double cr = cos(ang);
            double ci = sin(ang);
            
            sr += cr * psdph_r[k] - ci * psdph_i[k];
            si += cr * psdph_i[k] + ci * psdph_r[k];
        }
        
        // Coherence calculation remains in double
        double coh = sqrt(sr*sr + si*si) / sum_abs;
        
        if (coh > max_coh) {
            max_coh = coh;
            best_k_idx = t;
        }
    }
    
    *max_coh_out = max_coh;
    *best_k_idx_out = best_k_idx;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4) mexErrMsgIdAndTxt("Stamps:ps_topofit_mex:nrhs", "Inputs: ph, ph_patch, bperp_mat, n_trial_wraps");

    size_t n_ps = mxGetM(prhs[0]);
    size_t n_ifg = mxGetN(prhs[0]);
    int is_single = mxIsSingle(prhs[0]); 
    double *bperp_mat = mxGetPr(prhs[2]);
    double n_trial_wraps = mxGetScalar(prhs[3]);

    /* Outputs */
    plhs[0] = mxCreateDoubleMatrix(n_ps, 1, mxREAL); // K_ps
    plhs[1] = mxCreateDoubleMatrix(n_ps, 1, mxREAL); // C_ps
    plhs[2] = mxCreateDoubleMatrix(n_ps, 1, mxREAL); // coh_ps
    if (is_single) plhs[3] = mxCreateNumericMatrix(n_ps, n_ifg, mxSINGLE_CLASS, mxREAL);
    else plhs[3] = mxCreateDoubleMatrix(n_ps, n_ifg, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(n_ps, 1, mxREAL); // N_opt

    double *K_out = mxGetPr(plhs[0]);
    double *C_out = mxGetPr(plhs[1]);
    double *coh_out = mxGetPr(plhs[2]);
    void *ph_res_out = mxGetData(plhs[3]);
    double *N_opt_out = mxGetPr(plhs[4]);

    void *ph_data = mxGetData(prhs[0]);
    void *patch_data = mxGetData(prhs[1]);

    int limit = (int)ceil(8.0 * n_trial_wraps);
    int n_trials = 2 * limit + 1;
    double *trial_mult = (double*)calloc(n_trials, sizeof(double));
    for (int k = 0; k < n_trials; k++) trial_mult[k] = (double)(k - limit);

    double nan_val = mxGetNaN();

    #pragma omp parallel
    {
        /* Thread-local buffers */
        double *bperp_local = (double*)calloc(n_ifg, sizeof(double));
        double *psdph_r = (double*)calloc(n_ifg, sizeof(double));
        double *psdph_i = (double*)calloc(n_ifg, sizeof(double));
        double *weights = (double*)calloc(n_ifg, sizeof(double));
        int *valid_indices = (int*)calloc(n_ifg, sizeof(int));

        #pragma omp for
        for (size_t i = 0; i < n_ps; i++) {
            
            int n_valid = 0;
            double sum_abs = 0.0;
            double b_min = 1e30, b_max = -1e30;
            int has_invalid_ifg = 0; 

            /* 1. Data Extraction & STRICT Check */
            for (size_t j = 0; j < n_ifg; j++) {
                size_t idx = i + j * n_ps;
                double ph_r, ph_i, pt_r, pt_i;

                if (is_single) {
                    mxComplexSingle *d1 = (mxComplexSingle*)ph_data;
                    mxComplexSingle *d2 = (mxComplexSingle*)patch_data;
                    ph_r = d1[idx].real; ph_i = d1[idx].imag;
                    pt_r = d2[idx].real; pt_i = d2[idx].imag;
                } else {
                    mxComplexDouble *d1 = (mxComplexDouble*)ph_data;
                    mxComplexDouble *d2 = (mxComplexDouble*)patch_data;
                    ph_r = d1[idx].real; ph_i = d1[idx].imag;
                    pt_r = d2[idx].real; pt_i = d2[idx].imag;
                }

                // Initialize output residual to 0
                if (is_single) ((float*)ph_res_out)[idx] = 0.0f;
                else ((double*)ph_res_out)[idx] = 0.0;

                double rr = ph_r * pt_r + ph_i * pt_i;
                double ri = ph_i * pt_r - ph_r * pt_i;
                double mag = sqrt(rr*rr + ri*ri);

                /* STRICT CHECK: Matches MATLAB sum(psdph==0)==0 */
                if (mag <= 1e-9) {
                    has_invalid_ifg = 1;
                    break; // Fail immediately if ANY ifg is bad
                }

                psdph_r[n_valid] = rr;
                psdph_i[n_valid] = ri;
                weights[n_valid] = mag; 
                
                double b = bperp_mat[idx];
                bperp_local[n_valid] = b;
                
                if (b < b_min) b_min = b;
                if (b > b_max) b_max = b;
                
                sum_abs += mag;
                valid_indices[n_valid] = (int)j; 
                n_valid++;
            }

            /* If failed strict check, return NaN/0 exactly like MATLAB */
            if (has_invalid_ifg || n_valid == 0) {
                K_out[i] = nan_val;
                C_out[i] = 0.0; 
                coh_out[i] = 0.0; 
                N_opt_out[i] = 0.0;
                continue;
            }

            double bperp_range = b_max - b_min;
            if (bperp_range < 1e-6) bperp_range = 1.0; 

            /* 2. Coarse Search */
            double max_coh = -1.0;
            int best_k_idx = 0;
            double factor = (M_PI / 4.0) / bperp_range;

            // FORCED DOUBLE ACCUMULATION FOR STABILITY
            coarse_search_double_accum(n_trials, n_valid, trial_mult, bperp_local, factor, 
                                psdph_r, psdph_i, sum_abs, &max_coh, &best_k_idx);

            /* 3. Refinement (Weighted Least Squares) */
            double K0 = factor * trial_mult[best_k_idx];
            
            double sum_res_r = 0.0, sum_res_i = 0.0;
            for (int k = 0; k < n_valid; k++) {
                double ang = -(K0 * bperp_local[k]);
                double cr = cos(ang), ci = sin(ang);
                sum_res_r += psdph_r[k] * cr - psdph_i[k] * ci;
                sum_res_i += psdph_r[k] * ci + psdph_i[k] * cr;
            }
            double off_r = sum_res_r, off_i = -sum_res_i; 

            double numer = 0.0, denom = 0.0;
            for (int k = 0; k < n_valid; k++) {
                double ang = -(K0 * bperp_local[k]);
                double cr = cos(ang), ci = sin(ang);
                double rr = psdph_r[k] * cr - psdph_i[k] * ci;
                double ri = psdph_r[k] * ci + psdph_i[k] * cr;
                
                double yr = rr * off_r - ri * off_i;
                double yi = rr * off_i + ri * off_r;
                double ang_val = atan2(yi, yr);
                
                double wb = weights[k] * bperp_local[k];
                numer += wb * (weights[k] * ang_val);
                denom += wb * wb;
            }
            
            if (fabs(denom) > 1e-12) K0 += numer / denom;

            /* 4. Final Outputs */
            double fin_r = 0.0, fin_i = 0.0, fin_abs = 0.0;
            for (int k = 0; k < n_valid; k++) {
                double ang = -(K0 * bperp_local[k]);
                double cr = cos(ang), ci = sin(ang);
                double rr = psdph_r[k] * cr - psdph_i[k] * ci;
                double ri = psdph_r[k] * ci + psdph_i[k] * cr;
                fin_r += rr; fin_i += ri;
                fin_abs += sqrt(rr*rr + ri*ri);
                
                size_t out_idx = i + ((size_t)valid_indices[k]) * n_ps;
                double res_ang = atan2(ri, rr);
                
                // Final cast to single if necessary
                if (is_single) ((float*)ph_res_out)[out_idx] = (float)res_ang;
                else ((double*)ph_res_out)[out_idx] = res_ang;
            }

            K_out[i] = K0;
            C_out[i] = atan2(fin_i, fin_r);
            coh_out[i] = (fin_abs > 0) ? sqrt(fin_r*fin_r + fin_i*fin_i)/fin_abs : 0.0;
            N_opt_out[i] = (double)n_trials;
        }

        free(bperp_local); free(psdph_r); free(psdph_i);
        free(weights); free(valid_indices);
    } 

    free(trial_mult);
}