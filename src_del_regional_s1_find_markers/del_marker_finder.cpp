/*
    Given 
    
        a list of asPollinator-phased snp markers
        
    find the regions with large gap of certain size, e.g, >5 kb, to use them as "deletion" markers.
    
    2019-12-02 by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
*/
#include         <map>
#include      <string>
#include      <vector>
#include     <fstream>
#include     <sstream>
#include    <iostream>
#include   <algorithm>
#include    <string.h>
#include    <stdlib.h>
#include    <assert.h>
#include  <sys/stat.h>
#include    <dirent.h>
#include "split_string.h"
using namespace std;
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fmarker,
                 string* fsize,
                 int*    min_del_size,
                 string* outprefix,
                 bool*   verbose);
bool collect_contig_size(string                       sizefile, 
                         map<string, unsigned long>*  contigsize,
                         map<string, bool>*           contig_in_varlist);  
//
int main(int argc, char* argv[])
{
    if(argc < 9)
    {
        // g++ del_marker_finder.cpp split_string.cpp -O3 -o del_marker_finder
        cout << "\n   Function: find markers with large gaps in between, define that region as a \"deletion\" marker. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";        
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "             Note: \"deletion\" can be: real or artifacial resulted from the case that " << endl
             << "                   a homologous haplotig diverged from the other thus reads from the other unmapped."
             << endl;        
        cout << "\n   Usage: "
             << "del_marker_finder --marker snp_list.txt --min-del-size INT --chrsizes chrsizes.txt -o prefix_out" << endl << endl;
        return 1;
    }
    double startT= clock();
    // step 0. initialize from the cmd line
    cout << "Initialize input from cmd line: " << endl << endl;
    string fmarker   = "";
    string fsize     = "";    
    int min_del_size = 2000;
    string outprefix = "fun_";
    bool   verbose   = false; 
    if(!get_options(argc,
                    argv,
                    &fmarker,
                    &fsize,
                    &min_del_size,
                    &outprefix, 
                    &verbose) )
    {
        return 1;
    } 
    //
    vector<string> mkrfileinfo = split_string(fmarker, '/');
    std::stringstream ofilename;
    ofilename.str(""); 
    ofilename << "del_like_isize" << min_del_size << "_" << mkrfileinfo[mkrfileinfo.size()-1];
    ofstream ofp;
    ofp.open(ofilename.str().c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file " << ofilename.str() << endl;
        return 1;
    }
    //
    map<string, unsigned long> contigsize;
    map<string, bool> contig_in_varlist; // check if a contig has any variations    
    if(!collect_contig_size(fsize, &contigsize, &contig_in_varlist))
    {
        cout << "   Error: failed in collecting contig size info." << endl;
        return 1;
    }    
    // 
    ifstream ifp;
    ifp.open(fmarker.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << fmarker << endl;
        return 1;
    }
    string this_chr("");
    string last_chr(""); 
    unsigned long this_pos = 1;
    unsigned long last_pos = 1;
    bool last_chr_out = false;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        //proj_20190819_alt0	tig00000013_pilon	164	A	C	17.4923	3	1.0000	1
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5) 
        {
            cout << "   Warning: skipped insufficient info at line: " << line << endl;
            continue;
        }
        this_chr = lineinfo[1];
        this_pos = strtoul(lineinfo[2].c_str(), NULL, 0);
               
        if(last_chr.size()==0)
        {
            if(this_pos - 1 >= min_del_size)
            {
                // beginning of this_chr with last_chr null
                 ofp << this_chr  << "\t"
                     << 1         << "\t"
                     << this_pos  << endl; 
                 assert(contig_in_varlist.find(this_chr) != contig_in_varlist.end());
                 if(contig_in_varlist[this_chr] == false)
                 {
                     contig_in_varlist[this_chr] = true;
                 }                     
            }
        }else
        {
            if(last_chr.compare(this_chr) == 0)
            {
                // middle of current chr
                if(this_pos - last_pos >= min_del_size)
                {
                     ofp << this_chr  << "\t"
                         << last_pos  << "\t"
                         << this_pos  << endl; 
                     assert(contig_in_varlist.find(this_chr) != contig_in_varlist.end());
                     if(contig_in_varlist[this_chr] == false)
                     {
                         contig_in_varlist[this_chr] = true;
                     }                          
                }
            }else
            {
                // ending of last chr
                assert(contigsize.find(last_chr) != contigsize.end());
                if(contigsize[last_chr] - last_pos >= min_del_size)
                {
                     ofp << last_chr   << "\t"
                         << last_pos   << "\t"
                         << contigsize[last_chr]  << endl; 
                     assert(contig_in_varlist.find(this_chr) != contig_in_varlist.end());
                     if(contig_in_varlist[this_chr] == false)
                     {
                         contig_in_varlist[this_chr] = true;
                     }                          
                }
                /*
			last: proj_20190819_alt0  tig00000013_pilon   48630   A   C
			this: proj_20190819_alt0  tig00000026_pilon   13743   C   A          
                */
                // beginning of this_chr again with last_chr NOT null
                if(this_pos - 1 >= min_del_size)
                {
                     ofp << this_chr  << "\t"
                         << 1         << "\t"
                         << this_pos  << endl;       
                     assert(contig_in_varlist.find(this_chr) != contig_in_varlist.end());
                     if(contig_in_varlist[this_chr] == false)
                     {
                         contig_in_varlist[this_chr] = true;
                     }                                        
                }                              
            }
        }
        //
        last_chr = this_chr;
        last_pos = this_pos;        
    }
    // ending of last chr (reaching file end)
    if(contigsize[last_chr] - last_pos >= min_del_size)
    {
        ofp  << last_chr   << "\t"
             << last_pos   << "\t"
             << contigsize[last_chr]  << endl; 
        assert(contig_in_varlist.find(this_chr) != contig_in_varlist.end());
        if(contig_in_varlist[this_chr] == false)
        {
            contig_in_varlist[this_chr] = true;
        }              
    } 
    // also collect those in the assembly but not any variations defined.
    map<string, bool>::iterator citr;
    map<string, bool>::iterator citr_end;
    citr     = contig_in_varlist.begin();
    citr_end = contig_in_varlist.end();
    int iii  = 0;
    while(citr != citr_end)
    {
        if((*citr).second == false)
        {
            iii ++;
            cout << "   Info: " << (*citr).first << " without any variations collected as del-like. " << endl;
            ofp  << (*citr).first                << "\t"
                 << 1                            << "\t"
                 << contigsize[(*citr).first]    << endl;             
        }
        citr ++;
    }
    cout << "   Info: in total: " << iii << " such contigs. " << endl;
    //   
    ifp.close();
    ofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool collect_contig_size(string sizefile, map<string, unsigned long>* contigsize, map<string, bool>* contig_in_varlist)
{
    cout << "   Info: reading contig size info from "     << sizefile << endl;
    ifstream ifp;
    ifp.open(sizefile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open contig size file " << sizefile << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0) continue;
        if(line[0]=='#')     continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<2) continue;
        string        ctgid   = lineinfo[0];
        unsigned long ctgsize = strtoul(lineinfo[1].c_str(), NULL, 0);
        (*contigsize).insert(std::pair<string, unsigned long>(ctgid, ctgsize));
        (*contig_in_varlist).insert(std::pair<string, bool>(ctgid, false));
        //cout << ctgid << "\t" << ctgsize << endl;
    }
    if((*contigsize).size() == 0)
    {
        cout << "   Erorr: no info collected from contig size file. " << endl;
        return false;
    }else
    {
        cout << "   Info: " << (*contigsize).size() << " contig size info collected. " << endl;
    }
    ifp.close();
    //
    return true;
}
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fmarker,
                 string* fsize,
                 int*    min_del_size,
                 string* outprefix,
                 bool*   verbose)
{
    int  ic;
    // log cmdline
    cout << "   CMDline: ";
    ic = 0;
    while(ic < argc)
    {
        cout << argv[ic] << " ";
        ic ++;
    }
    cout << endl;
    // get values
    ic = 1;    
    while (ic < argc)
    {
        string optstr = (string)argv[ic];
        if(optstr.compare("--marker") == 0)
        {
            ic ++;
            *fmarker = (string)argv[ic];
            ifstream fp;
            fp.open((*fmarker).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: marker file provided: "         
                     << *fmarker << endl;
            }else
            {
                cout << "   Error: cannot open marker file "      
                     << *fmarker << endl;
                return false;
            }
        }else
        if(optstr.compare("--chrsizes") == 0)
        {
            ic ++;
            *fsize = (string)argv[ic];
            ifstream fp;
            fp.open((*fsize).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: chrsize file provided: "         
                     << *fsize << endl;
            }else
            {
                cout << "   Error: cannot open chrsize file "      
                     << *fsize << endl;
                return false;
            }
        }else        
        if(optstr.compare("--min-del-size") == 0)
        {
            ic ++;
            *min_del_size = atoi(argv[ic]);
            cout << "   Info: minimum deletion sizes provided as: "
                 << *min_del_size << endl;            
        }else
        if(optstr.compare("-o") == 0)
        {
            ic ++;
            *outprefix = (string)argv[ic];
            cout << "   Info: output files will be labeled with \""    
                 << *outprefix << "\"" << endl;
        }
        else              
        if(optstr.compare("-v") == 0)
        {
            *verbose = true;
            cout << "   Info: verbose mode turned on. "    
                 << endl;
        }
        else
        {
            cout << "   Warning: option "               << argv[ic] 
                 << " was not recognized and ignored."  << endl;
        }
        // next option
        ic ++;
    }
    // check necessary files - TODO
    return true;
}

