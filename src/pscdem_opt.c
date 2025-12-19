// *********************************************************************
// pscdem: Extract height for PS Candidates
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 10.03.2005
//
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 02/2010 AH Allow processing with no DEM
// 09/2015 AH Allow double precision
// 11/2025 Mingjia Li Optimize code and improve seek strategy
// ==============================================

#include <iostream>  
#include <fstream>
#include <vector>  
#include <string>  
#include <cmath>  
#include <cstdio>
#include <cstdlib>     
#include <cstring> 

using namespace std;

// Helper function to check file opening
void check_open(const ios& file, const char* name) {
    if (!file) {
        cerr << "pscdem: Error opening file " << name << endl;
        exit(1);
    }
}

int main(int argc, char *argv[]) {    

    // Optimization 1: Disable sync with stdio (Crucial for parsing large .ij files)
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);

    try {
        if (argc < 2) {	  
             cout << "Usage: pscdem pscands.1.ij pscands.1.hgt precision" << endl;
             throw "";
        }   
            
        const char *ijname;        
        if (argc < 3) ijname="pscands.1.ij";
        else ijname = argv[2];

        const char *outfilename;   
        if (argc < 4) outfilename="pscands.1.hgt";
        else outfilename = argv[3];

        const char *prec;
        if (argc < 5) prec="f";
        else prec = argv[4];

        ifstream parmfile (argv[1], ios::in);
        check_open(parmfile, argv[1]);

        ifstream psfile (ijname, ios::in); 
        cout << "opening " << ijname << "...\n";
        check_open(psfile, ijname);

        // Optimization 2: Large Output Buffer (1MB) to reduce write syscalls
        const size_t OUT_BUF_SIZE = 1024 * 1024; 
        char* out_buffer = new char[OUT_BUF_SIZE];
        ofstream outfile;
        outfile.rdbuf()->pubsetbuf(out_buffer, OUT_BUF_SIZE);
        outfile.open(outfilename, ios::out | ios::binary);
        check_open(outfile, outfilename);

        long width = 0;
        char ifgfilename[256];
        
        parmfile >> width;
        cout << "pscdem: width = " << width << "\n";	  
        
        string dummy; getline(parmfile, dummy);
        parmfile >> ifgfilename;

        // Optimization 3: Medium Input Buffer (64KB)
        // Helps capture "nearby" points without reading extremely large chunks.
        const size_t IN_BUF_SIZE = 64 * 1024;
        char* in_buffer = new char[IN_BUF_SIZE];
        ifstream ifgfile;
        ifgfile.rdbuf()->pubsetbuf(in_buffer, IN_BUF_SIZE);
        ifgfile.open(ifgfilename, ios::in | ios::binary);
        
        cout << "pscdem: opening " << ifgfilename << "...\n";

        int nodem_sw = 0;   
        if (!ifgfile.is_open()) {	    
            cout << "pscdem: Warning: cannot open " << ifgfilename << " - assumed all zeros" << endl; 
            nodem_sw = 1;
        } else {
            char header[32];
            long magic=0x59a66a95;
            ifgfile.read(header, 32);
            if (*reinterpret_cast<long*>(header) == magic)
                cout << "pscdem: sun raster file - skipping header\n";
            else 
                ifgfile.seekg(0, ios::beg); 
        }
        
        parmfile.close();

        long sizeofpixel; 
        if (prec[0]=='d') {
            sizeofpixel = sizeof(double);
            cout << "pscdem: input file specified as double precision\n";
        } else {
            sizeofpixel = sizeof(float);
            cout << "pscdem: input file specified as single precision\n";
        }

        // Buffer for a SINGLE pixel value (8 bytes is enough for double or float)
        char pixel_data[8] = {0}; 

        long pscid = 0;
        long x = 0;
        long y = 0;

        // Reverted Logic: Seek and Read
        while (psfile >> pscid >> y >> x) {
            
            if (nodem_sw == 0) {
                // Calculate absolute byte offset
                long long offset = (long long)(y * width + x) * sizeofpixel;
                
                ifgfile.seekg(offset, ios::beg);
                ifgfile.read(pixel_data, sizeofpixel);
            } else {
                // No DEM file: write zeros
                // (pixel_data is initialized to 0, but re-zeroing is safer inside loop logic)
                memset(pixel_data, 0, 8);
            }

            outfile.write(pixel_data, sizeofpixel);

            if (pscid % 100000 == 0 && pscid > 0)
                cout << "pscdem: " << pscid << " PS candidates processed\n";
        }   
        
        // Cleanup
        outfile.close();
        delete[] out_buffer;
        delete[] in_buffer;
        
    } catch(const char* str) {
        cout << str << "\n";
        return 999;
    } catch(...) {
        return 999;
    }

    return 0;
}