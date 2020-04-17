/* this function, given a number of cells with sequences of paternal and maternal patterns, 
   check number of potential swaps of alleles within each cells 
   
   Hequan Sun, MPIPZ, Email: sunhequan@gmail.com/sun@mpipz.mpg.de
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
#include                  <algorithm>
#include            "split_string.h"
//

//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << endl;
        cout << "   Given a number of cells with sequences of paternal and maternal patterns, output number of potential swap of haplotypes" << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "   (version 1.0 - compiled on " << buildString << ")"      << endl  << endl;
        cout << "   Usage: noise_checker list_of_pat_mat_sequences.file "   << endl;
        cout << endl;
        exit(1);
    }
    double startT= clock();
    //
    map<string, unsigned long> cellhapswaps;
    map<string, unsigned long> cellMkrCover;    
    //
    string filelistname = argv[1];
    ifstream ifp;
    ifp.open(filelistname.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << filelistname << endl;
        return 1;
    }
    while(ifp.good())
    {
        string line;
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]=='#') continue;
        
        ifstream pifp;
        pifp.open(line.c_str(), ios::in);
        if(!pifp.good())
        {
            cout << "   Warning: cannot open file " << line << endl;
            cout << "            skipped. "                 << endl;
            continue;
        }else
        {
            cout << "   Check: reading file " << line << endl;
        }
        while(pifp.good())
        {
            string line2;
            getline(pifp, line2);
            if(line2.size()==0) continue;
            if(line2[0]=='#') continue;            
            
            if( (line2.find("M")!=std::string::npos || line2.find("P")!=std::string::npos) &&
                line2.find("_u")==std::string::npos)
            {
                vector<string> line2info = split_string(line2, '\t');
                string pattern = line2info[0];
                string cellid  = line2info[1];
                //
                pattern.erase(std::remove(pattern.begin(), pattern.end(), 'u'), pattern.end()); // a sequence of "MMMPMMMP..."
                int markerCovered  = pattern.size();
                int occurrences    = 0;
                int            pos = 0;
                string      target = "PM";
                while ((pos = pattern.find(target, pos )) != std::string::npos) 
                {
                    ++ occurrences;
                    pos += target.length();
                }
                pos    = 0;
                target = "MP";
                while ((pos = pattern.find(target, pos )) != std::string::npos) 
                {
                    ++ occurrences;
                    pos += target.length();
                }
                //
                // if(occurrences > 0)
                cout << "   Check-cell: " <<  occurrences << " MP+PM found for " << markerCovered << " markers covered in cell " <<  cellid << endl;     
                // update counts
                map<string, unsigned long>::iterator tmpitr     = cellhapswaps.find(cellid);
                map<string, unsigned long>::iterator tmpitr_mkr = cellMkrCover.find(cellid);
                if(tmpitr == cellhapswaps.end() )
                {
                    cellhapswaps.insert(std::pair<string, unsigned long>(cellid, occurrences));
                    cellMkrCover.insert(std::pair<string, unsigned long>(cellid, markerCovered));
                }else
                {
                    (*tmpitr).second     += occurrences;
                    (*tmpitr_mkr).second += markerCovered;
                }
            }
            else
            if(line2.find("N")==std::string::npos && line2.find("_u")!=std::string::npos)
            {
                vector<string> line2info = split_string(line2, '\t');
                string pattern     = line2info[0];
                pattern.erase(std::remove(pattern.begin(), pattern.end(), 'u'), pattern.end()); // a sequence of "MMMPMMMP..."                
                int markerCovered  = pattern.size();
                string cellid      = line2info[1];                
                cout << "   Check-cell: 0 MP+PM found for " << markerCovered << " raw-markers covered in cell " <<  cellid << endl; 
            }
            else ;
        }
        pifp.close();
    }
    ifp.close();
    // output
    string outfile = "./haplotype_swaps_observed_in_cells.txt";
    ofstream ofp;
    ofp.open(outfile.c_str(), ios::out);
    map<string, unsigned long>::iterator tmpitr;
    map<string, unsigned long>::iterator tmpitr_end;
    tmpitr     = cellhapswaps.begin();
    tmpitr_end = cellhapswaps.end();
    while(tmpitr != tmpitr_end)
    {
        ofp << (*tmpitr).first << "\t" << (*tmpitr).second << "\t" << cellMkrCover[(*tmpitr).first] << endl;
        tmpitr ++; 
    }
    ofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return 0;
}
