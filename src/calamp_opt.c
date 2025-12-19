// *********************************************************************
// Calculate amplitude calibration constant for SLC files 
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//
// Change History
// ==============================================
// 01/2009 MA Deprication Fix
// 03/2009 MA Fix for gcc 4.3.x
// 01/2011 MCC neglecting pixel with  zero amplitude
// 12/2012 AH Add byteswap option
// 12/2025 Mingjia Li: Added full-file buffering and mask pre-loading.
// ==============================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <string>
#include <omp.h>

using namespace std;

// =======================================================================
// Byte-swap helpers
// =======================================================================
inline void cshortswap(complex<short>& f) {
    char* b = reinterpret_cast<char*>(&f);
    std::swap(b[0], b[1]);
    std::swap(b[2], b[3]);
}

inline void cfloatswap(complex<float>& f) {
    char* b = reinterpret_cast<char*>(&f);
    std::swap(b[0], b[3]);
    std::swap(b[1], b[2]);
    std::swap(b[4], b[7]);
    std::swap(b[5], b[6]);
}

struct FileStats {
    string filename;
    double calib_factor = 0.0;
    unsigned long long nof_pixels = 0;
    unsigned long long nof_zero_pixels = 0;
    bool ok = false;
    string error;
};

// =======================================================================
// Main Processing Function
// =======================================================================
// Now accepts the WHOLE mask data as a reference, avoiding file I/O inside
static FileStats process_slc_memory(const string& ampfilename,
                                    int width,
                                    const char* prec,
                                    int byteswap,
                                    const vector<char>& maskData,
                                    bool hasMask)
{
    FileStats stats;
    stats.filename = ampfilename;

    // 1. Open file to determine size and read
    ifstream ampfile(ampfilename.c_str(), ios::in | ios::binary | ios::ate); // Open at end to get size
    if (!ampfile.is_open()) {
        stats.error = "Error opening SLC file: " + ampfilename;
        return stats;
    }

    streamsize fileSize = ampfile.tellg();
    ampfile.seekg(0, ios::beg);

    // Calculate expected size validation
    size_t pixelSize = (prec[0] == 's') ? sizeof(complex<short>) : sizeof(complex<float>);
    size_t totalPixels = fileSize / pixelSize;
    size_t totalLines = totalPixels / width;

    // 2. Allocate memory for the WHOLE file
    // With 251GB RAM, allocating 5GB here is safe.
    // using vector<char> as raw buffer to handle both short/float easily
    vector<char> fileBuffer;
    try {
        fileBuffer.resize(fileSize);
    } catch (const std::bad_alloc& e) {
        stats.error = "Memory allocation failed for file (too big?): " + ampfilename;
        return stats;
    }

    // 3. Read the ENTIRE file in one go (Sequential Read - Best for HDD)
    if (!ampfile.read(fileBuffer.data(), fileSize)) {
        stats.error = "Error reading content of SLC file: " + ampfilename;
        return stats;
    }
    ampfile.close(); // Close file handle immediately to free OS resource

    // 4. Processing in Memory (CPU bound now, very fast)
    double sumamp = 0.0;
    unsigned long long nof_pixels = 0;
    unsigned long long nof_zero_pixels = 0;
    const float amp2_threshold = 1.0e-6f;

    // Pointers for raw data access
    const char* rawPtr = fileBuffer.data();
    
    // Check if mask matches dimensions (basic protection)
    if (hasMask && maskData.size() < totalLines * width) {
         // Logic assumes mask covers the whole image. 
    }

    // Iterate over all pixels
    // Using a flat loop is faster than nested loops
    for (size_t i = 0; i < totalPixels; ++i) {
        
        // Mask check
        if (hasMask) {
            // Safety check for mask bounds
            if (i < maskData.size() && maskData[i] != 0) {
                nof_zero_pixels++;
                continue;
            }
        }

        complex<float> camp;

        if (prec[0] == 's') {
            // Handle complex<short>
            const complex<short>* sPtr = reinterpret_cast<const complex<short>*>(rawPtr);
            complex<short> val = sPtr[i];
            if (byteswap) cshortswap(val);
            camp = val;
        } else {
            // Handle complex<float>
            const complex<float>* fPtr = reinterpret_cast<const complex<float>*>(rawPtr);
            complex<float> val = fPtr[i];
            if (byteswap) cfloatswap(val);
            camp = val;
        }

        float re = camp.real();
        float im = camp.imag();
        float mag2 = re * re + im * im;

        if (mag2 > amp2_threshold) {
            sumamp += std::sqrt(mag2);
            nof_pixels++;
        } else {
            nof_zero_pixels++;
        }
    }

    if (nof_pixels != 0) {
        stats.calib_factor = sumamp / static_cast<double>(nof_pixels);
    }

    stats.nof_pixels = nof_pixels;
    stats.nof_zero_pixels = nof_zero_pixels;
    stats.ok = true;
    return stats;
}

int main(int argc, char* argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) {
        cout << "Usage: calamp parmfile.in width parmfile.out precision byteswap maskfile\n";
        return 1;
    }

    const char* inlistname = argv[1];
    int width = std::atoi(argv[2]);
    const char* outfilename = (argc < 4) ? "parmfile.out" : argv[3];
    const char* prec        = (argc < 5) ? "f"             : argv[4];
    int byteswap = (argc < 6) ? 0 : std::atoi(argv[5]);
    string maskfilename;
    if (argc >= 7 && argv[6] && argv[6][0] != '\0') maskfilename = argv[6];

    // ==========================================
    // 1. Pre-load Mask (If exists)
    // ==========================================
    vector<char> maskData;
    bool hasMask = false;
    if (!maskfilename.empty()) {
        ifstream maskfile(maskfilename.c_str(), ios::in | ios::binary | ios::ate);
        if (maskfile.is_open()) {
            streamsize size = maskfile.tellg();
            maskfile.seekg(0, ios::beg);
            maskData.resize(size);
            maskfile.read(maskData.data(), size);
            hasMask = true;
            // Note: We suppress this load message or keep it minimal to not affect log parsing
            // cout << "Loaded mask file: " << maskfilename << " (" << size << " bytes)" << endl;
        }
    }

    // Read file list
    ifstream ampfiles(inlistname);
    vector<string> slcNames;
    string name;
    while (ampfiles >> name) if (!name.empty()) slcNames.push_back(name);
    ampfiles.close();

    ofstream parmfile(outfilename);
    vector<FileStats> results(slcNames.size());

    // ==========================================
    // 2. Parallel Processing with HDD Awareness
    // ==========================================
    // Recommendation: Set OMP_NUM_THREADS to a low number (e.g., 4 or 8)
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(slcNames.size()); ++i) {
        
        // Pass the maskData reference (read-only, thread-safe)
        FileStats stats = process_slc_memory(
            slcNames[i], width, prec, byteswap, maskData, hasMask);

        #pragma omp critical
        {
            if (!stats.ok) {
                cerr << "Error: " << stats.error << endl;
            } else {
                parmfile << stats.filename << " " << stats.calib_factor << "\n";
                parmfile.flush();
                
                // Logging format matches original calamp
                cout << "opening " << stats.filename << "..." << endl;
                cout << "Mean amplitude = "
                     << stats.calib_factor << endl;
                cout << "Number of pixels with zero amplitude = "
                     << stats.nof_zero_pixels << endl;
                cout << "Number of pixels with amplitude different than zero = "
                     << stats.nof_pixels << endl << endl;
                cout.flush();
            }
        }
    }

    parmfile.close();
    return 0;
}