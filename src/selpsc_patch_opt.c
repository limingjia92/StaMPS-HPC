// *********************************************************************
// Select PS Candidates
// Input are SLC's
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//  
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 07/2010 A0 Update Need for Speed on patches
// 08/2010 MA Code optimization
// 01/2010 MCC Drop low amplitudes
// 12/2012 AH Add byteswap and short options
// 12/2012 AH Correct mask processing
// 08/2017 AH add "else" in  mask processing
// 12/2025 Mingjia Li Optimize code (Block processing, I/O reduction)

// *********************************************************************
// Select PS Candidates (Optimized)
// Input are SLC's
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//  
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 07/2010 A0 Update Need for Speed on patches
// 08/2010 MA Code optimization
// 01/2010 MCC Drop low amplitudes
// 12/2012 AH Add byteswap and short options
// 12/2012 AH Correct mask processing
// 08/2017 AH add "else" in  mask processing
// 12/2025 Mingjia Li Optimize code (Block processing, I/O reduction, Log formatting)
// ==============================================

#include <iostream>
#include <fstream>
#include <string.h>
#include <complex>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cstdint>
#include <limits>
#include <iomanip>
#include <algorithm>

using namespace std;

// =======================================================================
// Swap functions
// =======================================================================
int cshortswap(complex<short>* f) {
    char* b = reinterpret_cast<char*>(f);
    complex<short> f2;
    char* b2 = reinterpret_cast<char*>(&f2);
    b2[0] = b[1];
    b2[1] = b[0];
    b2[2] = b[3];
    b2[3] = b[2];
    f[0] = f2;
    return 0;
}

int shortswap(short* f) {
    char* b = reinterpret_cast<char*>(f);
    short f2;
    char* b2 = reinterpret_cast<char*>(&f2);
    b2[0] = b[1];
    b2[1] = b[0];
    f[0] = f2;
    return 0;
}

int cfloatswap(complex<float>* f) {
    char* b = reinterpret_cast<char*>(f);
    complex<float> f2;
    char* b2 = reinterpret_cast<char*>(&f2);
    b2[0] = b[3];
    b2[1] = b[2];
    b2[2] = b[1];
    b2[3] = b[0];
    b2[4] = b[7];
    b2[5] = b[6];
    b2[6] = b[5];
    b2[7] = b[4];
    f[0] = f2;
    return 0;
}

int floatswap(float* f) {
    char* b = reinterpret_cast<char*>(f);
    float f2;
    char* b2 = reinterpret_cast<char*>(&f2);
    b2[0] = b[3];
    b2[1] = b[2];
    b2[2] = b[1];
    b2[3] = b[0];
    f[0] = f2;
    return 0;
}

int longswap(int32_t* f) {
    char* b = reinterpret_cast<char*>(f);
    int32_t f2;
    char* b2 = reinterpret_cast<char*>(&f2);
    b2[0] = b[3];
    b2[1] = b[2];
    b2[2] = b[1];
    b2[3] = b[0];
    f[0] = f2;
    return 0;
}

// =======================================================================
// Start of program 
// =======================================================================
int main(int argc, char* argv[]) {

    try {

        if (argc < 3) {
            cout << "Usage: selpsc parmfile patch.in pscands.1.ij pscands.1.da mean_amp.flt precision byteswap maskfile " << endl << endl;
            cout << "input parameters:" << endl;
            cout << "  parmfile (input) amplitude dispersion threshold" << endl;
            cout << "                   width of amplitude files (range bins)" << endl;
            cout << "                   SLC file names & calibration constants" << endl;
            cout << "  patch.in (input) location of patch in rg and az" << endl;
            cout << "  pscands.1.ij   (output) PS candidate locations" << endl;
            cout << "  pscands.1.da   (output) PS candidate amplitude dispersion" << endl << endl;
            cout << "  mean_amp.flt (output) mean amplitude of image" << endl << endl;
            cout << "  precision(input) s or f (default)" << endl;
            cout << "  byteswap   (input) 1 for to swap bytes, 0 otherwise (default)" << endl;
            cout << "  maskfile   (input)  mask rows and columns (optional)" << endl;
            throw "";
        }

        const char* ijname;
        if (argc < 4)
            ijname = "pscands.1.ij";
        else ijname = argv[3];

        char jiname[256]; // float format big endian for gamma
        strcpy(jiname, ijname);
        strcat(jiname, ".int");
        
        //MCC
        char ijname0[256]; //used to store PS with at least one amplitude =0.
        strcpy(ijname0, ijname);
        strcat(ijname0, "0");
        cout << "file name for zero amplitude PS: " << ijname0 << "\n";

        const char* daoutname;
        if (argc < 5)
            daoutname = "pscands.1.da";
        else daoutname = argv[4];

        const char* meanoutname;
        if (argc < 6)
            meanoutname = "mean_amp.flt";
        else meanoutname = argv[5];

        const char* prec;
        if (argc < 7)
            prec = "f";
        else prec = argv[6];

        int byteswap;
        if (argc < 8)
            byteswap = 0;
        else byteswap = atoi(argv[7]);

        const char* maskfilename;
        if (argc < 9)
            maskfilename = "";
        else maskfilename = argv[8];

        char referenceampfilename[256] = "0000";
        char referenceamp_exists = 0;
        if (argc < 10) {
        }
        else {
            ifstream referenceparmfile(argv[9], ios::in);
            if (!referenceparmfile.is_open()) {
                cout << "Error opening file " << argv[9] << "\n";
                throw "";
            }
            referenceparmfile >> referenceampfilename;
        }

        ifstream referenceampfile(referenceampfilename, ios::in);
        if (referenceampfile.is_open()) {
            referenceamp_exists = 1;
            cout << "opening " << referenceampfilename << "...\n";
        }

        ifstream parmfile(argv[1], ios::in);
        if (!parmfile.is_open()) {
            cout << "Error opening file " << argv[1] << "\n";
            throw "";
        }

        char line[256];
        int num_files = 0;

        int width = 0;
        float D_thresh = 0;
        int pick_higher = 0;

        parmfile >> D_thresh;
        cout << "dispersion threshold = " << D_thresh << "\n";
        float D_thresh_sq = D_thresh * D_thresh;
        if (D_thresh < 0) {
            pick_higher = 1;
        }

        parmfile >> width;
        cout << "width = " << width << "\n";
        parmfile.getline(line, 256);
        int savepos = parmfile.tellg();
        parmfile.getline(line, 256);
        while (!parmfile.eof()) {
            parmfile.getline(line, 256);
            num_files++;
        }
        
        parmfile.clear();
        parmfile.seekg(savepos);
        char ampfilename[256];
        ifstream* ampfile = new ifstream[num_files];
        float* calib_factor = new float[num_files];

        for (int i = 0; i < num_files; ++i) {
            parmfile >> ampfilename >> calib_factor[i];
            ampfile[i].open(ampfilename, ios::in | ios::binary);
            cout << "opening " << ampfilename << "...\n";

            if (!ampfile[i].is_open()) {
                cout << "Error opening file " << ampfilename << "\n";
                throw "";
            }

            char header[32];
            long magic = 0x59a66a95;
            ampfile[i].read(header, 32);
            if (*reinterpret_cast<long*>(header) == magic)
                cout << "sun raster file - skipping header\n";
            else ampfile[i].seekg(ios::beg);
        }

        // Optimization: Pre-calculate inverse calibration factors
        std::vector<float> calib_inv(num_files);
        for (int i = 0; i < num_files; ++i) {
            calib_inv[i] = (calib_factor[i] == 0.0f) ? 0.0f : (1.0f / calib_factor[i]);
        }

        parmfile.close();
        cout << "number of amplitude files = " << num_files << "\n";

        ifstream patchfile(argv[2], ios::in);
        if (!patchfile.is_open()) {
            cout << "Error opening file " << argv[2] << "\n";
            throw "";
        }

        int rg_start = 0;
        int rg_end = INT_MAX;
        int az_start = 0;
        int az_end = INT_MAX;
        patchfile >> rg_start;
        patchfile >> rg_end;
        patchfile >> az_start;
        patchfile >> az_end;
        patchfile.close();

        int patch_lines = az_end - az_start + 1;
        int patch_width = rg_end - rg_start + 1;

        const int sizeoffloat = 4;
        int sizeofelement;
        if (prec[0] == 's') {
            sizeofelement = sizeof(short);
        }
        else sizeofelement = sizeof(float);

        const int linebytes = width * sizeofelement * 2;
        const int patch_linebytes = patch_width * sizeofelement * 2;
        const int patch_amp_linebytes = patch_width * sizeofelement;

        filebuf* pbuf;
        long size;
        long numlines;

        pbuf = ampfile[0].rdbuf();
        size = pbuf->pubseekoff(0, ios::end, ios::in);
        pbuf->pubseekpos(0, ios::in);
        numlines = size / width / sizeofelement / 2;

        cout << "number of lines per file = " << numlines << "\n";
        cout << "patch lines = " << patch_lines << endl;
        cout << "patch width = " << patch_width << endl;

        ifstream maskfile(maskfilename, ios::in);
        char mask_exists = 0;
        if (maskfile.is_open()) {
            mask_exists = 1;
            cout << "opening " << maskfilename << "...\n";
        }

        ofstream ijfile(ijname, ios::out);
        ofstream jifile(jiname, ios::out);
        ofstream ijfile0(ijname0, ios::out);
        ofstream daoutfile(daoutname, ios::out);
        ofstream meanoutfile(meanoutname, ios::out);

        // Optimization: Adaptive buffer size
        const int LineInBufferMax = 128; // Increased buffer size
        int LineInBuffer = std::min(LineInBufferMax, patch_lines);

        // Memory allocation
        // Buffer stores [num_files] blocks. Each block has [LineInBufferMax] lines. Each line is [patch_linebytes] long.
        size_t bufferSize = (size_t)num_files * patch_linebytes * LineInBufferMax;
        char* buffer = new char[bufferSize];
        
        complex<float>* bufferf = reinterpret_cast<complex<float>*>(buffer);
        complex<short>* buffers = reinterpret_cast<complex<short>*>(buffer);

        char* maskbuffer = new char[patch_width * LineInBufferMax];
        // Initialize mask buffer (for safety, though it will be read over or ignored if no mask)
        memset(maskbuffer, 0, patch_width * LineInBufferMax);

        // Reference Amp Buffer
        char* referencebuffer = new char[patch_linebytes * LineInBufferMax];
        complex<float>* referencebufferf = reinterpret_cast<complex<float>*>(referencebuffer);
        complex<short>* referencebuffers = reinterpret_cast<complex<short>*>(referencebuffer);

        // Initialize reference buffer to default 1s
        for (int i = 0; i < patch_width * LineInBufferMax; i++) {
             if (prec[0] == 's') referencebuffers[i] = 1;
             else referencebufferf[i] = 1.0f;
        }

        long long pix_start_abs = (long long)(az_start - 1) * width + (rg_start - 1);
        long long pos_start = pix_start_abs * sizeofelement * 2;

        // Move all pointers to start position
        for (int i = 0; i < num_files; i++) {
            ampfile[i].seekg(pos_start, ios::beg);
        }
        if (mask_exists == 1) {
             maskfile.seekg(pix_start_abs, ios::beg);
        }
        if (referenceamp_exists == 1) {
            referenceampfile.seekg(pos_start, ios::beg);
        }

        int y_global = 0; // Line index within patch (0 to patch_lines-1)
        int pscid = 0;

        // =======================================================================
        // Optimized Processing Loop (Block Processing)
        // =======================================================================
        while (y_global < patch_lines) {
            
            // Determine how many lines to process in this block
            int lines_to_read = std::min(LineInBufferMax, patch_lines - y_global);
            
            // -------------------------------------------------------------------
            // 1. Read Phase: Load data into buffers
            // -------------------------------------------------------------------
            
            // Read Amplitudes
            for (int i = 0; i < num_files; i++) {
                char* file_ptr = &buffer[i * patch_linebytes * LineInBufferMax];
                for (int l = 0; l < lines_to_read; l++) {
                    ampfile[i].read(file_ptr + l * patch_linebytes, patch_linebytes);
                    if (l < lines_to_read) { 
                         ampfile[i].seekg(linebytes - patch_linebytes, ios::cur);
                    }
                }
            }

            // Read Mask
            if (mask_exists == 1) {
                for (int l = 0; l < lines_to_read; l++) {
                    maskfile.read(&maskbuffer[l * patch_width], patch_width);
                    maskfile.seekg(width - patch_width, ios::cur);
                }
            } else {
                 memset(maskbuffer, 0, patch_width * lines_to_read);
            }

            // Read Reference Amp
            if (referenceamp_exists == 1) {
                for (int l = 0; l < lines_to_read; l++) {
                    referenceampfile.read(&referencebuffer[l * patch_linebytes], patch_linebytes);
                    referenceampfile.seekg(linebytes - patch_linebytes, ios::cur);
                }
            } else {
                 for (int i = 0; i < patch_width * lines_to_read; i++) {
                     if (prec[0] == 's') referencebuffers[i] = 1;
                     else referencebufferf[i] = 1.0f;
                     
                     if (byteswap == 1) {
                         if (prec[0]=='s') cshortswap(&referencebuffers[i]); 
                         else cfloatswap(&referencebufferf[i]);
                     }
                }
            }

            // -------------------------------------------------------------------
            // 2. Compute Phase: Process the block
            // -------------------------------------------------------------------
            for (int l = 0; l < lines_to_read; l++) {
                int y = y_global + l; // Absolute line index in patch
                
                // Pointers for current line
                char* mask_line_ptr = &maskbuffer[l * patch_width];

                for (int x = 0; x < patch_width; x++) {
                    
                    float sumamp = 0;
                    float sumampsq = 0;
                    int amp_0 = 0;

                    // Fetch and Process Reference Amp for this pixel
                    complex<float> reference_amp_val;
                    int reference_idx = l * patch_width + x; 

                    if (prec[0] == 's') {
                        if (byteswap == 1 && referenceamp_exists == 1) { 
                            cshortswap(&referencebuffers[reference_idx]);
                        }
                        reference_amp_val = referencebuffers[reference_idx];
                    } else {
                        if (byteswap == 1 && referenceamp_exists == 1) {
                            cfloatswap(&referencebufferf[reference_idx]);
                        }
                        reference_amp_val = referencebufferf[reference_idx];
                    }

                    float abs_reference = abs(reference_amp_val);
                    if (abs_reference == 0) abs_reference = 1.0f;
                    float inv_reference = 1.0f / abs_reference;

                    // Loop over files
                    for (int i = 0; i < num_files; i++) {
                        complex<float> camp;
                        size_t idx = (size_t)i * patch_width * LineInBufferMax + (size_t)l * patch_width + x;

                        if (prec[0] == 's') {
                            if (byteswap == 1) {
                                cshortswap(&buffers[idx]);
                            }
                            camp = buffers[idx];
                        } else {
                            if (byteswap == 1) {
                                cfloatswap(&bufferf[idx]);
                            }
                            camp = bufferf[idx];
                        }

                        // Optimized math: use multiplication
                        float amp = abs(camp) * calib_inv[i] * inv_reference;

                        if (amp <= 0.00005f) {
                            amp_0 = 1;
                            sumamp = 0; 
                        } else {
                            sumamp += amp;
                            sumampsq += amp * amp;
                        }
                    }

                    meanoutfile.write(reinterpret_cast<char*>(&sumamp), sizeoffloat);

                    if (mask_line_ptr[x] == 0 && sumamp > 0) {
                        // Amplitude dispersion^2
                        // D_sq = num_files * sumampsq / (sumamp * sumamp) - 1;
                        double lhs_Dsq = (double)num_files * (double)sumampsq / ((double)sumamp * (double)sumamp) - 1.0;
                        
                        bool selected = false;
                        if (pick_higher == 0 && lhs_Dsq < D_thresh_sq) selected = true;
                        else if (pick_higher == 1 && lhs_Dsq >= D_thresh_sq) selected = true;

                        if (selected) {
                            if (amp_0 != 1) {
                                ++pscid;

                                ijfile << pscid << " " << (az_start - 1) + y << " " << (rg_start - 1) + x << "\n";
                                int32_t J = (rg_start - 1) + x;
                                int32_t I = (az_start - 1) + y;
                                longswap(&J);
                                longswap(&I);
                                jifile.write(reinterpret_cast<char*>(&J), sizeof(int32_t));
                                jifile.write(reinterpret_cast<char*>(&I), sizeof(int32_t));

                                float D_a = sqrt((float)lhs_Dsq);
                                daoutfile << D_a << "\n";
                            }
                            else {
                                ijfile0 << pscid << " " << (az_start - 1) + y << " " << (rg_start - 1) + x << "\n";
                            }
                        }
                    } // end mask check
                } // end x loop
            } // end processing loop (l)

            y_global += lines_to_read;

            // Progress reporting
            if ((y_global & 63) == 0 || y_global == patch_lines) {
                double pct = 100.0 * y_global / (double)patch_lines;
                cout << y_global << "/" << patch_lines
                     << " lines processed (" << fixed << setprecision(1) << pct << "%)"
                     << ", buffer contains " << lines_to_read << " lines"
                     << endl;
            }

        } // end while

        cout << endl << "Processing complete." << endl;

        ijfile.close();
        jifile.close();
        ijfile0.close();
        daoutfile.close();
        meanoutfile.close();
        
        if (mask_exists == 1) {
            maskfile.close();
        }
        if (referenceamp_exists == 1) {
            referenceampfile.close();
        }

        delete[] buffer;
        delete[] maskbuffer;
        delete[] referencebuffer;
        delete[] ampfile;
        delete[] calib_factor;
    }
    catch (char* str) {
        cout << str << "\n";
        return(999);
    }
    catch (...) {
        return(999);
    }

    return(0);

};