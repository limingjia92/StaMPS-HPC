/*
 * clap_filt_patch_mex.c
 * Batch Processor for PS Selection Filtering (Final Corrected Version)
 *
 * Revisions:
 * 1. Fixed compilation error: 'j_min_m' undeclared.
 * 2. Corrected fftshift/ifftshift logic for all window sizes.
 * 3. ADDED KERNEL NORMALIZATION to match MATLAB numerical behavior.
 * 4. Added explicit NaN handling for input grid.
 * 5. Optimized memory allocation.
 *
 * Compile: mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" clap_filt_patch_mex.c
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double r; double i; } Cpx;

/* Comparison for qsort */
int compare_doubles(const void *a, const void *b) {
    double arg1 = *(const double *)a;
    double arg2 = *(const double *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

/* 1D DFT using Lookup Table */
void dft_1d_lut(Cpx *in, Cpx *out, int n, Cpx *trig_tbl) {
    for (int k = 0; k < n; k++) {
        double sum_r = 0.0, sum_i = 0.0;
        int offset = k * n;
        for (int t = 0; t < n; t++) {
            double c = trig_tbl[offset + t].r;
            double s = trig_tbl[offset + t].i;
            sum_r += in[t].r * c - in[t].i * s;
            sum_i += in[t].r * s + in[t].i * c;
        }
        out[k].r = sum_r;
        out[k].i = sum_i;
    }
}

/* 2D DFT */
void dft_2d(Cpx *data, int n, Cpx *buf_vec, Cpx *buf_res, Cpx *tbl) {
    /* Rows */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) buf_vec[j] = data[i + j * n]; 
        dft_1d_lut(buf_vec, buf_res, n, tbl);
        for (int j = 0; j < n; j++) data[i + j * n] = buf_res[j];
    }
    /* Cols */
    for (int j = 0; j < n; j++) {
        size_t offset = j * n;
        for(int i=0; i<n; i++) buf_vec[i] = data[i + offset];
        dft_1d_lut(buf_vec, buf_res, n, tbl);
        for(int i=0; i<n; i++) data[i + offset] = buf_res[i];
    }
}

/* Separable Filter 2D 
 * Note: Performs Convolution (sum += img * kernel)
 * Matches MATLAB filter2/conv2 behavior
 */
void filter2_separable(double *img, int n, double *kernel_1d, int k_len, double *temp_buf) {
    int pad = k_len / 2;
    /* Columns (img -> temp_buf) */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < k_len; k++) {
                int r = i - pad + k;
                if (r >= 0 && r < n) sum += img[r + j * n] * kernel_1d[k];
            }
            temp_buf[i + j * n] = sum;
        }
    }
    /* Rows (temp_buf -> img) */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < k_len; k++) {
                int c = j - pad + k;
                if (c >= 0 && c < n) sum += temp_buf[i + c * n] * kernel_1d[k];
            }
            img[i + j * n] = sum;
        }
    }
}

/* Generic 2D Shift */
void fftshift_2d_generic(double *data, int n, int shift, double *tmp_buf) {
    memcpy(tmp_buf, data, n * n * sizeof(double));
    for(int j=0; j<n; j++) {
        int nj = (j + shift) % n;
        for(int i=0; i<n; i++) {
            int ni = (i + shift) % n;
            data[ni + nj*n] = tmp_buf[i + j*n];
        }
    }
}

void precompute_dft_table(int n, int direction, Cpx *tbl) {
    double theta_base = -direction * 2.0 * M_PI / (double)n;
    for (int k = 0; k < n; k++) {
        double theta_k = k * theta_base;
        for (int t = 0; t < n; t++) {
            double angle = theta_k * t;
            tbl[k * n + t].r = cos(angle);
            tbl[k * n + t].i = sin(angle);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 7) mexErrMsgIdAndTxt("clap_filt_patch_mex:args", "Inputs: ph_grid, grid_ij, alpha, beta, low_pass, osf, n_win, [B_1d]");

    /* 0. Input Parsing */
    const mwSize *dim_grid = mxGetDimensions(prhs[0]);
    int n_rows_grid = (int)dim_grid[0];
    int n_cols_grid = (int)dim_grid[1];
    int n_ifg = (mxGetNumberOfDimensions(prhs[0]) == 3) ? (int)dim_grid[2] : 1;
    mxComplexDouble *ph_grid = (mxComplexDouble*)mxGetData(prhs[0]);

    double *grid_ij = mxGetPr(prhs[1]); 
    size_t n_ps = mxGetM(prhs[1]);

    double alpha = mxGetScalar(prhs[2]);
    double beta = mxGetScalar(prhs[3]);
    double *low_pass = mxGetPr(prhs[4]); 
    
    int osf = (int)mxGetScalar(prhs[5]);
    int n_win = (int)mxGetScalar(prhs[6]);

    double *k_1d;
    int k_len = 0;
    if (nrhs >= 8) {
        k_1d = mxGetPr(prhs[7]);
        k_len = (int)mxGetM(prhs[7]) * (int)mxGetN(prhs[7]);
    } else {
        mexErrMsgIdAndTxt("clap_filt_patch_mex:args", "Please provide B_1d kernel as 8th input.");
    }

    /* Calculate Kernel Normalization Factor */
    /* This is critical to match MATLAB's numerical stability */
    double k_sum = 0.0;
    for(int k=0; k<k_len; k++) k_sum += k_1d[k];
    double k_norm_sq = k_sum * k_sum; // Separable filter means sum is squared

    /* Output */
    plhs[0] = mxCreateDoubleMatrix(n_ps, n_ifg, mxCOMPLEX);
    mxComplexDouble *ph_out = (mxComplexDouble*)mxGetData(plhs[0]);

    /* Precompute DFT Tables */
    Cpx *tbl_fwd = (Cpx*)malloc(n_win * n_win * sizeof(Cpx));
    Cpx *tbl_inv = (Cpx*)malloc(n_win * n_win * sizeof(Cpx));
    precompute_dft_table(n_win, 1, tbl_fwd);
    precompute_dft_table(n_win, -1, tbl_inv);

    /* Shift Amounts */
    int shift_fwd = n_win / 2;
    int shift_inv = (n_win + 1) / 2;

    #pragma omp parallel
    {
        /* Thread Local Buffers */
        int buf_sz = n_win * n_win;
        Cpx *patch_cpx = (Cpx*)malloc(buf_sz * sizeof(Cpx));
        double *patch_H = (double*)malloc(buf_sz * sizeof(double));
        double *patch_filt = (double*)malloc(buf_sz * sizeof(double));
        double *sort_buf = (double*)malloc(buf_sz * sizeof(double));
        double *shift_tmp = (double*)malloc(buf_sz * sizeof(double));
        Cpx *dft_vec = (Cpx*)malloc(n_win * sizeof(Cpx));
        Cpx *dft_res = (Cpx*)malloc(n_win * sizeof(Cpx));

        #pragma omp for
        for (size_t p = 0; p < n_ps; p++) {
            
            int r_idx = (int)grid_ij[p + 0*n_ps] - 1; 
            int c_idx = (int)grid_ij[p + 1*n_ps] - 1;

            // Row Bounds
            int i_min_m = (r_idx+1) - n_win/2;
            if (i_min_m < 1) i_min_m = 1;
            int i_max_m = i_min_m + n_win - 1;
            if (i_max_m > n_rows_grid) i_min_m = n_rows_grid - n_win + 1;
            
            // Col Bounds (FIXED: Added missing logic)
            int j_min_m = (c_idx+1) - n_win/2;
            if (j_min_m < 1) j_min_m = 1;
            int j_max_m = j_min_m + n_win - 1;
            if (j_max_m > n_cols_grid) j_min_m = n_cols_grid - n_win + 1;

            int r_start = i_min_m - 1;
            int c_start = j_min_m - 1;
            int p_rel_r = r_idx - r_start;
            int p_rel_c = c_idx - c_start;

            for (int k = 0; k < n_ifg; k++) {
                size_t page_offset = (size_t)k * (n_rows_grid * n_cols_grid);
                
                // 1. Extract Patch & Mask
                for (int c = 0; c < n_win; c++) {
                    for (int r = 0; r < n_win; r++) {
                        int src_r = r_start + r;
                        int src_c = c_start + c;
                        
                        // Mask logic
                        int dr = abs(r - p_rel_r);
                        int dc = abs(c - p_rel_c);
                        
                        if (dr < osf && dc < osf) {
                            patch_cpx[r + c*n_win].r = 0.0;
                            patch_cpx[r + c*n_win].i = 0.0;
                        } else {
                            if (src_r >= 0 && src_r < n_rows_grid && src_c >= 0 && src_c < n_cols_grid) {
                                size_t src_idx = page_offset + src_r + src_c * n_rows_grid;
                                double pr = ph_grid[src_idx].real;
                                double pi = ph_grid[src_idx].imag;
                                // Handle NaN input (Crucial for parity with MATLAB)
                                if (mxIsNaN(pr) || mxIsNaN(pi)) {
                                    patch_cpx[r + c*n_win].r = 0.0;
                                    patch_cpx[r + c*n_win].i = 0.0;
                                } else {
                                    patch_cpx[r + c*n_win].r = pr;
                                    patch_cpx[r + c*n_win].i = pi;
                                }
                            } else {
                                patch_cpx[r + c*n_win].r = 0.0;
                                patch_cpx[r + c*n_win].i = 0.0;
                            }
                        }
                    }
                }

                // 2. FFT
                dft_2d(patch_cpx, n_win, dft_vec, dft_res, tbl_fwd);

                // 3. Magnitude
                for(int i=0; i<buf_sz; i++) {
                    patch_H[i] = sqrt(patch_cpx[i].r*patch_cpx[i].r + patch_cpx[i].i*patch_cpx[i].i);
                }
                
                // 4. Smooth Response
                fftshift_2d_generic(patch_H, n_win, shift_fwd, shift_tmp);
                filter2_separable(patch_H, n_win, k_1d, k_len, patch_filt); // Output is in patch_H
                
                // *** NORMALIZATION STEP (CRITICAL) ***
                // Bring values back to ~1.0 scale to avoid precision issues in median/pow
                if (k_norm_sq > 1e-9) {
                    double k_inv = 1.0 / k_norm_sq;
                    for(int i=0; i<buf_sz; i++) patch_H[i] *= k_inv;
                }

                fftshift_2d_generic(patch_H, n_win, shift_inv, shift_tmp);

                // 5. Median Normalization
                memcpy(sort_buf, patch_H, buf_sz * sizeof(double));
                qsort(sort_buf, buf_sz, sizeof(double), compare_doubles);
                double meanH;
                if (buf_sz % 2 == 0) meanH = 0.5 * (sort_buf[buf_sz/2 - 1] + sort_buf[buf_sz/2]);
                else meanH = sort_buf[buf_sz/2];

                // 6. Apply Filter G
                double invH = (meanH > 1e-12) ? 1.0/meanH : 0.0;
                for (int i=0; i<buf_sz; i++) {
                    if (meanH > 1e-12) patch_H[i] *= invH;
                    
                    double val = pow(patch_H[i], alpha) - 1.0;
                    if (val < 0) val = 0.0;
                    
                    // Low pass is already aligned (DC at corners)
                    double G = val * beta + low_pass[i]; 
                    
                    patch_cpx[i].r *= G;
                    patch_cpx[i].i *= G;
                }

                // 7. IFFT
                dft_2d(patch_cpx, n_win, dft_vec, dft_res, tbl_inv);

                // 8. Store Result
                double scale = 1.0 / (double)buf_sz;
                size_t res_idx = p_rel_r + p_rel_c * n_win;
                
                ph_out[p + k * n_ps].real = patch_cpx[res_idx].r * scale;
                ph_out[p + k * n_ps].imag = patch_cpx[res_idx].i * scale;
            }
        }

        free(patch_cpx); free(patch_H); free(patch_filt); free(sort_buf);
        free(shift_tmp); free(dft_vec); free(dft_res);
    }
    free(tbl_fwd); free(tbl_inv);
}