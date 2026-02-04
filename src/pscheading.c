// *********************************************************************
// pscheading: Calculate Heading and Incidence angles for PS Candidates
// ---------------------------------------------------------------------
// AUTHOR    : Mingjia Li
// ---------------------------------------------------------------------
// WRITTEN   : Jan 2026

#include <iostream>  
#include <fstream>  
#include <vector>
#include <string>  
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ==========================================
// Helper: Byte Swapping (Big Endian <-> Little Endian)
// ==========================================
void swap_float(float *f) {
    char *ptr = reinterpret_cast<char*>(f);
    char temp;
    temp = ptr[0]; ptr[0] = ptr[3]; ptr[3] = temp;
    temp = ptr[1]; ptr[1] = ptr[2]; ptr[2] = temp;
}

// ==========================================
// Structs & Utils
// ==========================================
struct DemParams {
    long width = 0;
    long nlines = 0;
    double corner_lat = 0.0;
    double corner_lon = 0.0;
    double post_lat = 0.0;
    double post_lon = 0.0;
};

void check_open(const ios& file, const char* name) {
    if (!file) {
        cerr << "Error: Cannot open file " << name << endl;
        exit(1);
    }
}

DemParams read_dem_params(const char* filename) {
    ifstream file(filename);
    check_open(file, filename);

    DemParams params;
    string line, key;
    
    // Parse loop compatible with GAMMA par files
    while (file >> key) {
        if (key == "width:") file >> params.width;
        else if (key == "nlines:") file >> params.nlines;
        else if (key == "corner_lat:") file >> params.corner_lat;
        else if (key == "corner_lon:") file >> params.corner_lon;
        else if (key == "post_lat:") file >> params.post_lat;
        else if (key == "post_lon:") file >> params.post_lon;
        getline(file, line); 
    }
    file.close();
    return params;
}

// ==========================================
// Main
// ==========================================
int main(int argc, char *argv[]) {    

    // Optimization: Disable sync with stdio
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);

    try {
        // 1. Argument Parsing
        // Expected: ./pscheading pscheading.in pscands.1.ll pscands.1.head pscands.1.inc
        if (argc != 5) {	  
            cerr << "Usage: " << argv[0] << " pscheading.in pscands.1.ll pscands.1.head pscands.1.inc" << endl;
            return 1;
        }

        const char *parmfile_name = argv[1];
        const char *ll_name       = argv[2];
        const char *head_out_name = argv[3];
        const char *inc_out_name  = argv[4];

        // 2. Read Paths from pscheading.in
        ifstream pfile(parmfile_name);
        check_open(pfile, parmfile_name);
        
        string dem_par_path, theta_path, phi_path;
        pfile >> dem_par_path >> theta_path >> phi_path;
        pfile.close();

        // 3. Read DEM Parameters
        DemParams dem = read_dem_params(dem_par_path.c_str());
        
        // Basic Validation
        if (dem.width <= 0 || dem.post_lon == 0 || dem.post_lat == 0) {
            cerr << "Error: Invalid DEM parameters (Width=" << dem.width << ")" << endl;
            return 1;
        }

        // 4. Open Binary Files (Input)
        ifstream f_theta(theta_path.c_str(), ios::in | ios::binary);
        check_open(f_theta, theta_path.c_str());

        ifstream f_phi(phi_path.c_str(), ios::in | ios::binary);
        check_open(f_phi, phi_path.c_str());

        ifstream f_ll(ll_name, ios::in | ios::binary);
        check_open(f_ll, ll_name);

        // Optimization: Set Input Buffers (64KB)
        const size_t IN_BUF_SIZE = 64 * 1024; 
        vector<char> buf_theta(IN_BUF_SIZE);
        vector<char> buf_phi(IN_BUF_SIZE);
        f_theta.rdbuf()->pubsetbuf(buf_theta.data(), IN_BUF_SIZE);
        f_phi.rdbuf()->pubsetbuf(buf_phi.data(), IN_BUF_SIZE);

        // 5. Open Binary Files (Output)
        ofstream f_out_head(head_out_name, ios::out | ios::binary);
        check_open(f_out_head, head_out_name);

        ofstream f_out_inc(inc_out_name, ios::out | ios::binary);
        check_open(f_out_inc, inc_out_name);

        // Optimization: Set Output Buffers (1MB)
        const size_t OUT_BUF_SIZE = 1024 * 1024; 
        char* buf_head = new char[OUT_BUF_SIZE];
        char* buf_inc  = new char[OUT_BUF_SIZE];
        f_out_head.rdbuf()->pubsetbuf(buf_head, OUT_BUF_SIZE);
        f_out_inc.rdbuf()->pubsetbuf(buf_inc, OUT_BUF_SIZE);

        // 6. Processing Loop
        float lon, lat;
        float theta_val, phi_val;
        float head_res, inc_res;
        long x, y;
        long long seek_addr;
        long count = 0;
        long out_of_bounds = 0;

        // Loop reads 4 bytes (lon) then 4 bytes (lat)
        while (f_ll.read(reinterpret_cast<char*>(&lon), sizeof(float))) {
            if (!f_ll.read(reinterpret_cast<char*>(&lat), sizeof(float))) break;

            // SWAP INPUT: GAMMA pscands.1.ll is Big Endian
            swap_float(&lon);
            swap_float(&lat);

            // NEAREST NEIGHBOR LOOKUP
            // Calculate pixel index based on corner and spacing
            // Use round() to handle non-aligned coordinates (nearest pixel)
            x = (long)round((lon - dem.corner_lon) / dem.post_lon);
            y = (long)round((lat - dem.corner_lat) / dem.post_lat);

            // Boundary Checks
            if (x >= 0 && x < dem.width && y >= 0 && y < dem.nlines) {
                
                seek_addr = (long long)(y * dem.width + x) * sizeof(float);

                // Read Theta (Big Endian)
                f_theta.seekg(seek_addr, ios::beg);
                f_theta.read(reinterpret_cast<char*>(&theta_val), sizeof(float));
                swap_float(&theta_val); // Convert to Host(LE) for calculation

                // Read Phi (Big Endian)
                f_phi.seekg(seek_addr, ios::beg);
                f_phi.read(reinterpret_cast<char*>(&phi_val), sizeof(float));
                swap_float(&phi_val); // Convert to Host(LE) for calculation

                // Calculations
                // Heading = -phi/pi*180 - 180
                head_res = -phi_val / (float)M_PI * 180.0f - 180.0f;
                
                // Incidence = 90 - theta/pi*180
                inc_res = 90.0f - theta_val / (float)M_PI * 180.0f;

            } else {
                // Handle points outside DEM coverage (set to 0 or NaN)
                out_of_bounds++;
                head_res = 0.0f;
                inc_res = 0.0f;
            }

            // SWAP OUTPUT: Write as Big Endian (Standard for GAMMA/StaMPS chain)
            swap_float(&head_res);
            swap_float(&inc_res);

            f_out_head.write(reinterpret_cast<char*>(&head_res), sizeof(float));
            f_out_inc.write(reinterpret_cast<char*>(&inc_res), sizeof(float));

            count++;
            if (count % 1000000 == 0 && count > 0)
                cout << "pscheading: " << count << " PS candidates processed\n";
        }

        cout << "Done. Total: " << count << ", Out of bounds: " << out_of_bounds << endl;

        // Cleanup
        delete[] buf_head;
        delete[] buf_inc;

    } catch(const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 999;
    }

    return 0;
}