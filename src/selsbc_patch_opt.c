// *********************************************************************
// Select SB Candidates
// Input are SLC's
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 04.08.2003
//  
// Change History
// ==============================================
// 03/2009 MA  Fix for gcc 4.3.x
// 12/2012 LIG Speed up processing by reading blocks
// 12/2012 DB  Merge with developers version (SVN)
// 02/2018 DB  Read 10 lines at the time to handle larger temporal datasets
// 08/2018 GJF Prevent integer overflow for buffer size
// 12/2025 Mingjia Li Optimize code (Block processing, I/O reduction, Log formatting)
// ==============================================

#include <iostream>  
using namespace std;
     
#include <fstream>  
using namespace std;

#include <string.h>  
using namespace std;
     
#include <complex>  
using namespace std;
     
#include <vector>  
using namespace std;
     
#include <cmath>  
using namespace std;
     
#include <cstdio>
using namespace std;     

#include <cstdlib>     
using namespace std;     

#include <climits>     
using namespace std;     

#include <cstdint>     
using namespace std;     

#include <limits>
#include <cstring>  
#include <iomanip>

// =======================================================================
// Swap functions
// =======================================================================
int cshortswap( complex<short>* f )
{
  char* b = reinterpret_cast<char*>(f);
  complex<short> f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[1];
  b2[1] = b[0];
  b2[2] = b[3];
  b2[3] = b[2];
  f[0]=f2;
  return 0;
}

int cfloatswap( complex<float>* f )
{
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
  f[0]=f2;
  return 0;
}

int longswap( int32_t* f )
{
  char* b = reinterpret_cast<char*>(f);
  int32_t f2;
  char* b2 = reinterpret_cast<char*>(&f2);
  b2[0] = b[3];
  b2[1] = b[2];
  b2[2] = b[1];
  b2[3] = b[0];
  f[0]=f2;
  return 0;
}

// =======================================================================
// Start of program 
// =======================================================================
int main(int  argc, char *argv[] ) {   

try {
 
  if (argc < 3)
  {	  
     cout << "Usage: selsbc parmfile patch.in pscands.1.ij pscands.1.da mean_amp.flt precision byteswap maskfile" << endl << endl;
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
     
  const char *ijname;          
  if (argc < 4) 
     ijname="pscands.1.ij";
  else ijname = argv[3];   

  char jiname[256]; // float format big endian for gamma
  strcpy (jiname,ijname);
  strcat (jiname,".int");
     
  const char *daoutname;       
  if (argc < 5) 
     daoutname="pscands.1.da";
  else daoutname = argv[4];   
     
  const char *meanoutname;
  if (argc < 6) 
     meanoutname="mean_amp.flt";
  else meanoutname = argv[5];   
  
  (void)meanoutname;

  const char *prec;
  if (argc < 7)
     prec="f";
  else prec = argv[6];

  int byteswap;
  if (argc < 8)
     byteswap=0;
  else byteswap = atoi(argv[7]);
  
  const char *maskfilename;   
  if (argc < 9) 
     maskfilename="";
  else maskfilename = argv[8];   
     
  ifstream parmfile (argv[1], ios::in);
  if (! parmfile.is_open()) 
  {	  
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
  float D_thresh_sq = D_thresh*D_thresh;
  if (D_thresh<0) { 
      pick_higher=1;
  }

  parmfile >> width;
  cout << "width = " << width << "\n";	  
  parmfile.getline(line,256);
  int savepos=parmfile.tellg();
  parmfile.getline(line,256);
  while (! parmfile.eof())
  {
      parmfile.getline(line,256);
      num_files++;
  }    

  parmfile.clear();
  parmfile.seekg(savepos);
  char ampfilename[256];
  ifstream* ampfile   = new ifstream[num_files];
  float* calib_factor = new float[num_files];
      
  for (int i=0; i<num_files; ++i) 
  {
    parmfile >> ampfilename >> calib_factor[i];
    ampfile[i].open (ampfilename, ios::in|ios::binary);
    cout << "opening " << ampfilename << " [file " << i << "]...\n";

    if (! ampfile[i].is_open())
    {	    
        cout << "Error opening file " << ampfilename << "\n"; 
	throw "";
    }

    char header[32];
    long magic=0x59a66a95;
    ampfile[i].read(header,32);
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

  ifstream patchfile (argv[2], ios::in);
  if (! patchfile.is_open()) 
  {	  
      cout << "Error opening file " << argv[2] << "\n"; 
      throw "";
  }    

  int rg_start=0;
  int rg_end=INT_MAX;
  int az_start=0;
  int az_end=INT_MAX;
  patchfile >> rg_start;
  patchfile >> rg_end;
  patchfile >> az_start;
  patchfile >> az_end;
  patchfile.close();

  int sizeofelement; 
  if (prec[0]=='s')
  {
      sizeofelement = sizeof(short);
  }else sizeofelement = sizeof(float);

  const size_t linebytes = width*sizeofelement*2;  // bytes per line in amplitude files (SLCs)


  filebuf *pbuf;
  long size;
  long numlines;

  // get pointer to associated buffer object
  pbuf=ampfile[0].rdbuf();

  // get file size using buffer's members
  size=pbuf->pubseekoff (0,ios::end,ios::in);
  pbuf->pubseekpos (0,ios::in);
  numlines=size/width/sizeof(float)/2;

  cout << "number of lines per file = " << numlines << "\n";	  
  
  ifstream maskfile (maskfilename, ios::in);
  char mask_exists = 0;
  if (maskfile.is_open()) 
  {	 
      cout << "opening " << maskfilename << "...\n";
      mask_exists=1;
  }    
  
  ofstream ijfile(ijname,ios::out);
  ofstream jifile(jiname,ios::out);
  ofstream daoutfile(daoutname,ios::out);
 
  // Optimization: Adaptive buffer size
  const int LineInBufferMax = 128;     
  int LineInBuffer = std::min(LineInBufferMax, az_end - (az_start-1));  // Initial block size

  size_t BufferSize;

  cout << num_files << " files, " << linebytes << " line bytes, " << LineInBuffer << " lines in the buffer\n";

  BufferSize = num_files*linebytes*LineInBufferMax;
  char* buffer = new char[BufferSize];
    
  complex<float>* bufferf = reinterpret_cast<complex<float>*>(buffer);
  complex<short>* buffers = reinterpret_cast<complex<short>*>(buffer);

  char* maskline = new char[width*LineInBufferMax];

  for (int x=0; x<width*LineInBuffer; x++) // for each pixel in range
  {
      maskline[x] = 0;
  }

  int y=az_start-1;                             // Current azimuth line (0-based)
  int pscid=0;                                  // PS candidate ID number

  long long pix_start; 				// start position in azimuth in values
  long long pos_start; 				// start position of azimuth in bytes
  pix_start= (long long)(y)*(long long)(width);
  pos_start = (long long)(y) * (long long)(linebytes);
  
    
  for ( int i=0; i<num_files; i++)              // read in first line from each amp file
  {
      ampfile[i].seekg (pos_start, ios::beg);
      ampfile[i].read (reinterpret_cast<char*>(&buffer[i*linebytes*LineInBuffer]), linebytes*LineInBuffer);
  } 
     
  if (mask_exists==1) 
  {
     maskfile.seekg (pix_start, ios::beg);      // set pointer to start of patch in mask file
  }

  while (! ampfile[1].eof() && y < az_end) 
  {
     if (mask_exists==1) 
     {
         maskfile.read (maskline, width*LineInBuffer);
     }    
     if (y >= az_start-1)
     {
      for(int j=0; j<LineInBuffer; j++) //for each line in azimuth within buffer
      {
       for (int x=rg_start-1; x<rg_end; x++) // for each pixel in range
       {
        
        // Check external mask if it exists
        if (mask_exists==1 && maskline[j*width+x] != 0) {
            continue;
        }
        
        float sumamp = 0;
        float sumampdiffsq = 0;
        int i;

        for (i=0; i<num_files/2; i++)        // for each amp file
	    {

           complex<float> camp1, camp2; // get amp value
           if (prec[0]=='s')
           {
               if (byteswap == 1)
               {
                  cshortswap(&buffers[(i*2)*width*LineInBuffer+j*width+x]);
                  cshortswap(&buffers[(i*2+1)*width*LineInBuffer+j*width+x]);
               }
               camp1=buffers[(i*2)*width*LineInBuffer+j*width+x];
               camp2=buffers[(i*2+1)*width*LineInBuffer+j*width+x];
           }
           else
           {
               camp1=bufferf[(i*2)*width*LineInBuffer+j*width+x]; // get amp value
               camp2=bufferf[(i*2+1)*width*LineInBuffer+j*width+x]; // get amp value
               if (byteswap == 1)
               {
                  cfloatswap(&camp1);
                  cfloatswap(&camp2);
               }
           }

           // Optimization: Use pre-calculated inverse calibration factors
           float amp1 = abs(camp1) * calib_inv[i*2];
           float amp2 = abs(camp2) * calib_inv[i*2+1];

           sumamp+=amp1;
           sumamp+=amp2;
           sumampdiffsq+=(amp1-amp2)*(amp1-amp2);
        }

        if (maskline[j*width+x]==0 && sumamp > 0)
        {
            // Optimization: Use double precision for threshold check to match previous logic
            double nPairs = (double)(num_files/2);
            double nImgs  = (double)num_files;
            double lhs = (double)sumampdiffsq / nPairs;           
            double rhs = (double)D_thresh_sq * ( (double)sumamp / nImgs ) * ( (double)sumamp / nImgs );

            if ( (pick_higher==0 && lhs <  rhs) ||
                 (pick_higher==1 && lhs >= rhs) ) {

                double meanA  = (double)sumamp / nImgs;
                double rmsD   = std::sqrt( (double)sumampdiffsq / nPairs );
                double Da     = (meanA>0.0) ? (rmsD / meanA) : numeric_limits<double>::quiet_NaN();

                ++pscid;
                ijfile << pscid << " " << y << " " << x << "\n";

                int32_t J=x, I=y;
                longswap(&J); longswap(&I);
                jifile.write(reinterpret_cast<char*>(&J), sizeof(int32_t));
                jifile.write(reinterpret_cast<char*>(&I), sizeof(int32_t));

                daoutfile << Da << "\n";
            }
        }

       } // x++
               y++;
           }//j++
     } // endif y>az_start-1
      
     // Update buffer size for next read
     LineInBuffer = (az_end - y) > LineInBufferMax ? LineInBufferMax : (az_end - y);

     for ( int i=0; i<num_files; i++)           // read in next line from each amp file
     {
        ampfile[i].read (reinterpret_cast<char*>(&buffer[i*linebytes*LineInBuffer]), linebytes*LineInBuffer);
     } 
     
     // Progress reporting (optimized to reduce I/O spam)
     long processed = y - (az_start - 1);              
     long total     = az_end - (az_start - 1);         

     if ( (processed & 63) == 0 ) {                    
        double pct = 100.0 * processed / (double)total;
        cout << processed << "/" << total
             << " lines processed (" << fixed << setprecision(1) << pct << "%)"
             << ", buffer contains " << LineInBuffer << " lines"
             << std::endl;                            
     }
  }  //while
  
  ijfile.close();
  daoutfile.close();
  jifile.close();
  
  if (mask_exists==1) 
  {	  
      maskfile.close();
  }    
  
  }


  catch( char * str ) {	
     cout << str << "\n";
     return(999);
  }   

  return(0);
       
};