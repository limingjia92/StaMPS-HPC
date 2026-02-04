// *********************************************************************
// psclonlat: Extract lon/lat for PS Candidates
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 10.03.2005
//
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 11/2025 Mingjia Li Optimize code and improve seek strategy
// ==============================================


#include <iostream>  
#include <fstream>  
#include <vector>
#include <string>  
#include <cstdlib>     

using namespace std;

// Helper to check file opening
void check_open(const ios& file, const char* name) {
    if (!file) {
        cerr << "Error opening file " << name << endl;
        exit(1);
    }
}

int main(int argc, char *argv[]) {    

    // Optimization 1: Disable sync with stdio (Speeds up parsing large .ij files)
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);

    try {
        if (argc < 2) {	  
            cout << "Usage: psclonlat parmfile pscands.1.ij pscands.1.ll" << endl;
            throw "";
        }   
            
        const char *ijname;          
        if (argc < 3) ijname="pscands.1.ij";
        else ijname = argv[2];

        const char *outfilename;     
        if (argc < 4) outfilename="pscands.1.ll";
        else outfilename = argv[3];

        ifstream parmfile(argv[1], ios::in);
        check_open(parmfile, argv[1]);

        ifstream psfile(ijname, ios::in); 
        cout << "opening " << ijname << "...\n";
        check_open(psfile, ijname);

        // Optimization 2: 1MB Buffer for Output (Reduces write syscalls)
        const size_t OUT_BUF_SIZE = 1024 * 1024; 
        char* out_buffer = new char[OUT_BUF_SIZE];
        ofstream outfile;
        outfile.rdbuf()->pubsetbuf(out_buffer, OUT_BUF_SIZE);
        outfile.open(outfilename, ios::out | ios::binary);
        check_open(outfile, outfilename);

        int num_files = 2;
        long width = 0;
        char ifgfilename[256];
        
        parmfile >> width;
        cout << "width = " << width << "\n";	  
        string dummy; getline(parmfile, dummy); 

        vector<ifstream> ifgfiles(num_files);
        
        // Optimization 3: Set larger internal read buffers for input files
        // This allows the OS/Stream to cache "nearby" data without forcing a full row read.
        // 64KB is a sweet spot for random seeks on HDD.
        const size_t IN_BUF_SIZE = 64 * 1024; 
        vector<vector<char>> in_bufs(num_files, vector<char>(IN_BUF_SIZE));

        for (int i=0; i<num_files; ++i) {
            parmfile >> ifgfilename;
            
            // Set buffer BEFORE opening
            ifgfiles[i].rdbuf()->pubsetbuf(in_bufs[i].data(), IN_BUF_SIZE);
            ifgfiles[i].open(ifgfilename, ios::in | ios::binary);
            
            cout << "opening " << ifgfilename << "...\n";
            check_open(ifgfiles[i], ifgfilename);

            char header[32];
            long magic=0x59a66a95;
            ifgfiles[i].read(header, 32);
            if (*reinterpret_cast<long*>(header) == magic)
                cout << "sun raster file - skipping header\n";
            else 
                ifgfiles[i].seekg(0, ios::beg); 
        }
        parmfile.close();

        long pscid = 0;
        long x = 0;
        long y = 0;
        float val;

        // Reverted Loop: Seek and Read
        // This is more efficient when points are sparse (Width >> Points_per_row)
        while (psfile >> pscid >> y >> x) {
            
            long long xyaddr = (long long)(y * width + x) * sizeof(float);

            for (int i = 0; i < num_files; ++i) {
                ifgfiles[i].seekg(xyaddr, ios::beg);
                ifgfiles[i].read(reinterpret_cast<char*>(&val), sizeof(float));
                outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
            }

            if (pscid % 1000000 == 0 && pscid > 0)
                cout << "psclonlat: " << pscid << " PS candidates processed\n";

        }

        outfile.close();
        delete[] out_buffer;
        
    } catch(const char* str) {
        cout << str << "\n";
        return 999;
    } catch(...) {
        return 999;
    }

    return 0;
}