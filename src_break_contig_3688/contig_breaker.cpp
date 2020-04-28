/* this function breaks contigs of very likely chimeric.

   Hequan Sun, MPIPZ Email:sunhequan@gmail.com/sun@mpipz.mpg.de
*/
#include                   <stdio.h>
#include                  <string.h>
#include                    <string>
#include                       <map>
#include                  <stdlib.h>
#include                  <iostream>
#include                   <sstream>
#include                   <fstream>
#include                    <time.h>
#include                  <assert.h>
#include                   <iomanip>
#include                <sys/stat.h>
#include                  <dirent.h>
#include                 <algorithm>
#include            "split_string.h"
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << endl;
        cout << "   Given a marker/consensus file, break contig 3688 to 2 parts. "   << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "   (version 1.0 - compiled on " << buildString << ")"       << endl  << endl;
        cout << "   Usage: contig_breaker file INT "                         << endl;
        cout << "          if INT=0, consensus file so contig-id at first column; else marker file so contig-id at second column" << endl;
        cout << endl;
        exit(1);
    }
    double startT= clock();
    /*
      consensus format:
      tig00003688_pilon	176254	G	1	0	0	1	0	0	0	1

      marker format:
      proj_20190819_alt0	tig00003688_pilon	33943	C	G	222	141	0.5382	1
    */
    int concol = atoi(argv[2]);
    if(concol != 0) concol = 1;

    // contigs to be separated
    map<string, unsigned long> sepcontigs;
    string contig             = "tig00003688_pilon";
    unsigned long separatepos = 2050000;
    sepcontigs.insert(std::pair<string, unsigned long>(contig, separatepos) );
    contig="tig00003826_pilon";
    separatepos=220000; // 210579
    sepcontigs.insert(std::pair<string, unsigned long>(contig, separatepos) );

    string infile = (string)argv[1];
    fstream ifp;
    ifp.open(infile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << infile << endl;
        return 1;
    }
    vector<string> infileinfo = split_string(infile, '.');
    string ofilename = infileinfo[0] + "_update.txt";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file " << ofilename << endl;
        return 1;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0) continue;
        if(line[0]=='#') continue;

        vector<string> lineinfo = split_string(line, '\t');
        contig                  = lineinfo[concol];

        if(sepcontigs.find(contig) != sepcontigs.end() && contig.find("pilons")==std::string::npos)
        {
            separatepos = sepcontigs[contig];
            unsigned long  pos = strtoul(lineinfo[concol+1].c_str(), NULL, 0);
            if(pos < separatepos)
            {
                for(int ii=0; ii< lineinfo.size(); ii++)
                {
                    if(ii == concol)
                    {
                       ofp << lineinfo[ii] << "sA\t";
                    }else
                    {
                       if(ii < lineinfo.size()-1)
                       {
                           ofp << lineinfo[ii] << "\t";
                       }else
                       {
                           ofp << lineinfo[ii] << endl;
                       }
                    }
                }
            }else
            {
                for(int ii=0; ii< lineinfo.size(); ii++)
                {
                    if(ii == concol)
                    {
                       ofp << lineinfo[ii] << "sB\t";
                    }else
                    {
                       if(ii < lineinfo.size()-1)
                       {
                           ofp << lineinfo[ii] << "\t";
                       }else
                       {
                           ofp << lineinfo[ii] << endl;
                       }
                    }
                }
            }
        }else
        {
            ofp << line << endl;
        }

    }
    ifp.close();
    ofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return 0;
}
