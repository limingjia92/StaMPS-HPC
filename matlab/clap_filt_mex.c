/*
 * clap_filt_mex.c (Highly Optimized Version)
 * * Optimizations:
 * 1. Separable Convolution (Splits 2D 7x7 kernel into two 1D 7-tap filters).
 * 2. Trigonometric Lookup Tables (Avoids sin/cos calls in DFT).
 * 3. Pre-allocated memory buffers to reduce heap contention.
 *
 * Compile: mex -R2018a CFLAGS="$CFLAGS -fopenmp -O3 " LDFLAGS="$LDFLAGS -fopenmp" clap_filt_mex.c
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
    /* trig_tbl layout: [k*t] where index = k*n + t */
    for (int k = 0; k < n; k++) {
        double sum_r = 0.0;
        double sum_i = 0.0;
        int offset = k * n;
        
        for (int t = 0; t < n; t++) {
            /* exp(-j*theta) = cos - j*sin */
            /* Table stores: r=cos, i=-sin (for forward) or i=sin (for inverse) */
            double c = trig_tbl[offset + t].r;
            double s = trig_tbl[offset + t].i;
            
            sum_r += in[t].r * c - in[t].i * s;
            sum_i += in[t].r * s + in[t].i * c;
        }
        out[k].r = sum_r;
        out[k].i = sum_i;
    }
}

/* 2D DFT using LUT */
void dft_2d(Cpx *data, int n, Cpx *buf_in, Cpx *buf_out, Cpx *trig_tbl) {
    /* Rows */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) buf_in[j] = data[i + j * n]; 
        dft_1d_lut(buf_in, buf_out, n, trig_tbl);
        for (int j = 0; j < n; j++) data[i + j * n] = buf_out[j];
    }
    
    /* Cols */
    for (int j = 0; j < n; j++) {
        size_t offset = j * n;
        for (int i = 0; i < n; i++) buf_in[i] = data[i + offset];
        dft_1d_lut(buf_in, buf_out, n, trig_tbl);
        for (int i = 0; i < n; i++) data[i + offset] = buf_out[i];
    }
}

/* Standard fftshift */
void fftshift_2d(double *data, int n) {
    int mid = n / 2;
    for (int j = 0; j < mid; j++) {
        for (int i = 0; i < mid; i++) {
            size_t idx1 = i + j * n;
            size_t idx4 = (i + mid) + (j + mid) * n;
            size_t idx2 = (i + mid) + j * n;
            size_t idx3 = i + (j + mid) * n;
            double tmp = data[idx1]; data[idx1] = data[idx4]; data[idx4] = tmp;
            tmp = data[idx2]; data[idx2] = data[idx3]; data[idx3] = tmp;
        }
    }
}

/* * Separable Filter 2D (Convolution) 
 * Kernel B is assumed to be symmetric Gaussian (7x1) derived from gausswin(7).
 * This replaces the O(K^2) per pixel with O(2*K).
 */
void filter2_separable(double *img, int n, double *kernel_1d, int k_len, double *temp_buf) {
    int pad = k_len / 2;
    
    /* 1. Convolve Columns (Vertical) -> Store in temp_buf */
    /* Accessing column-wise is cache friendly for 'img' which is col-major */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < k_len; k++) {
                int r = i - pad + k;
                if (r >= 0 && r < n) {
                    sum += img[r + j * n] * kernel_1d[k];
                }
            }
            temp_buf[i + j * n] = sum;
        }
    }

    /* 2. Convolve Rows (Horizontal) on temp_buf -> Store back in img */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < k_len; k++) {
                int c = j - pad + k;
                if (c >= 0 && c < n) {
                    /* Row i is contiguous in memory? No, Col j is.
                       temp_buf[i + c*n] accesses widely separated memory.
                       But for small N (32), it fits in L1 cache mostly.
                    */
                    sum += temp_buf[i + c * n] * kernel_1d[k];
                }
            }
            img[i + j * n] = sum;
        }
    }
}

/* Helper to precompute DFT tables */
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
    if (nrhs != 8) mexErrMsgIdAndTxt("clap_filt_mex:args", "8 inputs required.");

    /* Inputs */
    size_t n_i = mxGetM(prhs[0]);
    size_t n_j = mxGetN(prhs[0]);
    mxComplexDouble *ph = (mxComplexDouble*)mxGetData(prhs[0]);
    
    double alpha = mxGetScalar(prhs[1]);
    double beta  = mxGetScalar(prhs[2]);
    int n_win    = (int)mxGetScalar(prhs[3]); // e.g. 24
    int n_pad    = (int)mxGetScalar(prhs[4]); // e.g. 8
    double *low_pass = mxGetPr(prhs[5]);
    
    /* Input 6 is the 2D Kernel B. We need to extract 1D kernel from it.
       The MATLAB code constructs B = gausswin(7)*gausswin(7)'.
       We can take the first column of B (assuming it's symmetric/rank-1) and normalize?
       Wait, simpler: Let's assume the user PASSES the 1D gausswin vector as input 6 
       OR we extract the middle column. 
       Since we can't easily change the MATLAB caller in this snippet, let's extract the middle column/row.
       B is k_rows x k_cols. 
    */
    double *B_2d = mxGetPr(prhs[6]); 
    int k_size = (int)mxGetM(prhs[6]); // Assumed square 7
    
    /* Extract 1D Kernel (Middle Column) and Sqrt it? 
       No, if B = v * v', then v = sqrt(diag(B)) or similar.
       However, technically B_2d is passed. 
       Let's create a 1D kernel. Since B = G * G', sum(B) = sum(G)*sum(G').
       Actually, taking the middle column of B gives G * G_middle. 
       Let's just compute the marginal sum or take the sqrt of the diagonal if passed properly.
       
       SAFE HACK: The kernel is small (7). We can just take the middle column and normalize it 
       such that sum(kernel_1d)^2 = sum(B_2d). 
       
       Better yet: User B_2d is separable. 
       Let's just take the middle column B(:, ceil(k/2)).
       If B = v*v', then col_mid = v * v_mid. 
       We need 'v'. 
       v = col_mid / sqrt(B(mid, mid)).
    */
    
    double *kernel_1d = (double*)malloc(k_size * sizeof(double));
    int mid_idx = k_size / 2;
    double center_val = B_2d[mid_idx + mid_idx * k_size]; // B(mid, mid)
    if (center_val > 0) {
        double scale = 1.0 / sqrt(center_val);
        for(int k=0; k<k_size; k++) {
            kernel_1d[k] = B_2d[k + mid_idx * k_size] * scale;
        }
    } else {
        // Fallback (Identity)
        memset(kernel_1d, 0, k_size*sizeof(double));
        kernel_1d[mid_idx] = 1.0;
    }

    double *wind_func_tpl = mxGetPr(prhs[7]);

    int n_win_ex = n_win + n_pad; // 32
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(n_i, n_j, mxCOMPLEX);
    mxComplexDouble *ph_out = (mxComplexDouble*)mxGetData(plhs[0]);

    int n_inc = n_win / 4; 
    int n_win_i = (int)ceil((double)n_i / n_inc) - 3;
    int n_win_j = (int)ceil((double)n_j / n_inc) - 3;

    /* Precompute Trig Tables for DFT (Forward and Inverse) */
    Cpx *tbl_fwd = (Cpx*)malloc(n_win_ex * n_win_ex * sizeof(Cpx));
    Cpx *tbl_inv = (Cpx*)malloc(n_win_ex * n_win_ex * sizeof(Cpx));
    precompute_dft_table(n_win_ex, 1, tbl_fwd);
    precompute_dft_table(n_win_ex, -1, tbl_inv);
    
    /* OpenMP Parallel Loop */
    #pragma omp parallel
    {
        /* Thread Local Allocations */
        size_t buf_size = n_win_ex * n_win_ex;
        Cpx *ph_bit = (Cpx*)calloc(buf_size, sizeof(Cpx));
        double *H = (double*)calloc(buf_size, sizeof(double));
        double *H_sort = (double*)calloc(buf_size, sizeof(double));
        
        /* Scratch buffers */
        Cpx *dft_buf_in = (Cpx*)calloc(n_win_ex, sizeof(Cpx));
        Cpx *dft_buf_out = (Cpx*)calloc(n_win_ex, sizeof(Cpx));
        double *filt_buf = (double*)calloc(buf_size, sizeof(double));

        #pragma omp for collapse(1) schedule(dynamic)
        for (int ix1 = 0; ix1 < n_win_i; ix1++) {
            int i1 = ix1 * n_inc; 
            int i2 = i1 + n_win - 1;
            int i_shift = 0;
            if (i2 >= n_i) { i_shift = i2 - n_i + 1; i2 = n_i - 1; i1 = n_i - n_win; }

            for (int ix2 = 0; ix2 < n_win_j; ix2++) {
                int j1 = ix2 * n_inc;
                int j2 = j1 + n_win - 1;
                int j_shift = 0;
                if (j2 >= n_j) { j_shift = j2 - n_j + 1; j2 = n_j - 1; j1 = n_j - n_win; }

                /* 1. Extract Patch */
                memset(ph_bit, 0, buf_size * sizeof(Cpx));
                for (int c = 0; c < n_win; c++) {
                    for (int r = 0; r < n_win; r++) {
                        size_t src_idx = (i1 + r) + (j1 + c) * n_i;
                        size_t dst_idx = r + c * n_win_ex; 
                        ph_bit[dst_idx].r = ph[src_idx].real;
                        ph_bit[dst_idx].i = ph[src_idx].imag;
                        // NaN check omitted for speed if data is clean, or keep if needed
                        if (isnan(ph_bit[dst_idx].r)) { ph_bit[dst_idx].r=0; ph_bit[dst_idx].i=0; }
                    }
                }

                /* 2. FFT2 (LUT) */
                dft_2d(ph_bit, n_win_ex, dft_buf_in, dft_buf_out, tbl_fwd);
                
                /* 3. Magnitude */
                for (int k = 0; k < buf_size; k++) {
                    H[k] = sqrt(ph_bit[k].r*ph_bit[k].r + ph_bit[k].i*ph_bit[k].i);
                }
                
                /* 4. Filter (Separable + Shift) */
                fftshift_2d(H, n_win_ex);
                // Use Separable convolution here:
                filter2_separable(H, n_win_ex, kernel_1d, k_size, filt_buf);
                fftshift_2d(H, n_win_ex);
                
                /* 5. Median Normalize (Partial sort could be better, but Qsort is robust) */
                memcpy(H_sort, H, buf_size * sizeof(double));
                qsort(H_sort, buf_size, sizeof(double), compare_doubles);
                double meanH;
                if (buf_size % 2 == 0) meanH = 0.5 * (H_sort[buf_size/2 - 1] + H_sort[buf_size/2]);
                else meanH = H_sort[buf_size/2];

                if (meanH > 1e-9) {
                    double invH = 1.0 / meanH;
                    for (int k = 0; k < buf_size; k++) H[k] *= invH;
                }
                
                /* 6. Alpha/Beta/LowPass */
                for (int k = 0; k < buf_size; k++) {
                    double val = pow(H[k], alpha) - 1.0;
                    if (val < 0) val = 0.0;
                    double G = val * beta + low_pass[k];
                    ph_bit[k].r *= G;
                    ph_bit[k].i *= G;
                }
                
                /* 7. IFFT2 (LUT) */
                dft_2d(ph_bit, n_win_ex, dft_buf_in, dft_buf_out, tbl_inv);
                
                /* 8. Accumulate */
                double norm_factor = 1.0 / ((double)n_win_ex * (double)n_win_ex);
                for (int c = 0; c < n_win; c++) {
                    if (c < j_shift) continue; 
                    for (int r = 0; r < n_win; r++) {
                        if (r < i_shift) continue;
                        
                        size_t tpl_idx = (r - i_shift) + (c - j_shift) * n_win;
                        double win_val = wind_func_tpl[tpl_idx];
                        
                        if (win_val > 0.0) {
                            size_t res_idx = r + c * n_win_ex;
                            size_t out_idx = (i1 + r) + (j1 + c) * n_i;
                            
                            double val_r = ph_bit[res_idx].r * norm_factor * win_val;
                            double val_i = ph_bit[res_idx].i * norm_factor * win_val;
                            
                            #pragma omp atomic
                            ph_out[out_idx].real += val_r;
                            #pragma omp atomic
                            ph_out[out_idx].imag += val_i;
                        }
                    }
                }
            } 
        } 
        
        free(ph_bit); free(H); free(H_sort);
        free(dft_buf_in); free(dft_buf_out); free(filt_buf);
    }
    
    free(tbl_fwd); free(tbl_inv); free(kernel_1d);
}