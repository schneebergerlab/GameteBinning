/*
    Given 
    
        a list of del_marker_finder defined bed del-like regions, and
        the samtools depth defined depth,
        
    define del-like regions as hap, hom or repeat.
    
    2019-12-28 by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
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
                 string* fregion,
                 string* fdepth,
                 int*    hap_hom_cutoff,
                 string* outprefix,
                 bool*   verbose);
bool get_del_bed(string  fregion,
                 map<string, unsigned long>*  regions);
bool get_del_value(string fdepth, 
                 map<string, unsigned long >* regions);
bool output_del_depth(map<string, unsigned long>   regions, 
                 string outprefix);                                  
int hap_hom_cutoff = 200;                 
//
int main(int argc, char* argv[])
{
    if(argc < 9)
    {
        // g++ del_depth_finder.cpp split_string.cpp -O3 -o del_depth_finder
        cout << "\n   Function: find depth for del-like regions; define them as hap, hom (including repeat-like). ";
        const char *buildString = __DATE__ ", " __TIME__ ".";        
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "             Note: \"deletion\" can be: real or artifacial resulted from the case that " << endl
             << "                   a homologous haplotig diverged from the other thus reads from the other unmapped."
             << endl;        
        cout << "\n   Usage: "
             << "del_depth_finder --region del_like.bed --depth del_like_depth.txt --hap-hom-cutoff INT -o prefix_out" 
             << endl << endl;
        cout << "          Note: --hap-hom-cutoff determines which depths would be defined as hap/hom region. \n"
             << "          You should determine this value with histogram of depth in genome sequencing. \n"
             << "          For example, cutoff = 1/2 * genome-wide average sequencing depth. " << endl << endl;
        return 1;
    }
    double startT= clock();
    // step 0. initialize from the cmd line
    cout << "Initialize input from cmd line: " << endl << endl;
    string fregion   = "";
    string fdepth    = "";    
    string outprefix = "fun_";
    bool   verbose   = false; 
    if(!get_options(argc,
                    argv,
                    &fregion,
                    &fdepth,
                    &hap_hom_cutoff,
                    &outprefix, 
                    &verbose) )
    {
        return 1;
    } 
    // get del-like regions
    map<string, unsigned long> regions;
    if(!get_del_bed(fregion, &regions) )
    {
        cout << "   Error: cannot read del-like beds. " << endl;
        return 1;
    }
    // get average depth in del-like regions
    if(!get_del_value(fdepth, &regions))
    {
        cout << "   Error: cannot read del-like depth. " << endl;
        return 1;
    }
    // output 
    if(!output_del_depth(regions, outprefix))
    {
        cout << "   Error: cannot write output. " << endl;
        return 1;
    }    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool output_del_depth(map<string, unsigned long> regions, string outprefix)
{
    string ofilename = outprefix + "_del_like_interval_avg_depth.bed";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: failed in writing output to file " << ofilename << endl;
        return false;
    }
    // ofp << "#aheader:contig\tstart\tend\tdepth\ttype_[0, 200]x:hap_high,(200, 250]x:hap_low, (250, 300]x:hom_low,(300, inf)x:hom_high" << endl;
    ofp << "#aheader:contig\tstart\tend\tdepth\ttype_" 
        << "["    << "0"                               << ", " << hap_hom_cutoff                    << "]x:hap_high;"
        << "("    << hap_hom_cutoff                    << ", " << hap_hom_cutoff + hap_hom_cutoff/4 << "]x:hap_low;"
        << "("    << hap_hom_cutoff + hap_hom_cutoff/4 << ", " << hap_hom_cutoff + hap_hom_cutoff/2 << "]x:hom_low;"
        << "("    << hap_hom_cutoff + hap_hom_cutoff/2 << ", " << "inf"                             << "]x:hom_high"
        << endl;    
    map<string, unsigned long>::iterator oitr;
    map<string, unsigned long>::iterator oitr_end;
    oitr     = regions.begin();
    oitr_end = regions.end();
    while(oitr != oitr_end)
    {
        ofp << (*oitr).first << "\t" << (*oitr).second << "\t";
        if((*oitr).second <= hap_hom_cutoff)
        {
            ofp << "hap_high" << endl;
        }else
        if((*oitr).second > hap_hom_cutoff && (*oitr).second <= hap_hom_cutoff+hap_hom_cutoff/4)
        {
            ofp << "hap_low"  << endl;
        }else        
        if((*oitr).second > hap_hom_cutoff+hap_hom_cutoff/4 && (*oitr).second <= hap_hom_cutoff+hap_hom_cutoff/2)
        {
            ofp << "hom_low"  << endl;
        }else
        if((*oitr).second > hap_hom_cutoff+hap_hom_cutoff/2)
        {
            ofp << "hom_high" << endl;
        }
        oitr ++;
    }
    ofp.close();
    //
    std::stringstream ss;
    ss.str("");
    ss << "sort -k1,1 -k2,2n " << ofilename << " > " << outprefix << "_del_like_interval_avg_depth_sorted.bed";
    system(ss.str().c_str());
    //cout << "   Info: you can sort the file with " << endl << endl
    //     << "         " << ss.str() << endl;
    //
    return true;
}
//
bool get_del_value(string fdepth, map<string, unsigned long >* regions)
{
    ifstream ifp;
    ifp.open(fdepth.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: failed in reading " << fdepth << endl;
        return false;
    }
    string               key = "";
    unsigned long region_dep = 0;
    unsigned long region_cnt = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        
        // tig00000013_pilon   5118    239
        vector<string> lineinfo = split_string(line, '\t');
        string        ctg = lineinfo[0];
        unsigned long pos = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long dep = strtoul(lineinfo[2].c_str(), NULL, 0);
        std::stringstream ss;
        ss.str("");        
        if(key.size()==0)
        {
            ss << pos - 1;            
            key         = ctg + "\t" + ss.str() + "\t";
            region_dep += dep;
            region_cnt ++;
            // cout << "   Check: key initialized as " << key << endl;
        }else
        {
            region_dep += dep;    
            region_cnt ++;        
            ss << pos;
            string tmpkey = key + ss.str();
            if((*regions).find(tmpkey) != (*regions).end())
            {
                if(region_cnt == 0) region_cnt = 1; // caution 0
                (*regions)[tmpkey] = region_dep / region_cnt;
                cout << "   Found: " << tmpkey << " with depth " << region_dep / region_cnt << " from " << region_cnt << " lines" << endl;
                // next region
                key = "";
                region_dep = 0;
                region_cnt = 0;
            }else
            {
                // cout << "   Check: tmpkey \"" << tmpkey << "\" not found in collected regions. "<< endl;
            }
        }
    }
    return true;
}
//
bool get_del_bed(string fregion, map<string, unsigned long>* regions)
{
    ifstream ifp;
    ifp.open(fregion.c_str(), ios::in);
    if(!ifp.good())
    {
       cout << "   Error: failed in reading " << fregion << endl;
       return false;
    }
    int regionnum = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        if(line.size() < 3) 
        {
            cout << "   Warning: skipped insufficient line info at " << line << endl;
            continue;
        }
        
        string        ctg = lineinfo[0]; 
        string        sta = lineinfo[1]; 
        string        end = lineinfo[2]; 
        string        key = ctg + "\t" + sta + "\t" + end;
        regionnum ++;     
        
        map<string, unsigned long >::iterator ctgitr;
        ctgitr = (*regions).find(key);
        (*regions).insert(std::pair<string,unsigned long>(key, 0));
        // tig00000026_pilon	37715	105094
    }
    cout << "   Info: " << (*regions).size() << " regions collected." << endl;
    ifp.close();
    return true;
}
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fregion,
                 string* fdepth,
                 int*    hap_hom_cutoff,
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
        if(optstr.compare("--region") == 0)
        {
            ic ++;
            *fregion = (string)argv[ic];
            ifstream fp;
            fp.open((*fregion).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: region file provided: "         
                     << *fregion << endl;
            }else
            {
                cout << "   Error: cannot open region file "      
                     << *fregion << endl;
                return false;
            }
        }else
        if(optstr.compare("--depth") == 0)
        {
            ic ++;
            *fdepth = (string)argv[ic];
            ifstream fp;
            fp.open((*fdepth).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: depth file provided: "         
                     << *fdepth << endl;
            }else
            {
                cout << "   Error: cannot open depth file "      
                     << *fdepth << endl;
                return false;
            }
        }else
        if(optstr.compare("--hap-hom-cutoff") == 0)
        {
            ic ++;
            *hap_hom_cutoff = atoi(argv[ic]);
            cout << "   Info: hap-hom depth cutoff provided: " << *hap_hom_cutoff << endl;
            if(*hap_hom_cutoff < 0)
            {
                cout << "   Error: hap-hom depth cutoff must be > 0. " << endl;
                return false;
            }
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

