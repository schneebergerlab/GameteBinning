/* this function, given genetic map file, and initial fasta file, make contigs into chrs.
   Written by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
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
bool reverse_complement(string* kmer, string * rckmer);
bool read_genetic_map(string gmfile, vector<string>* gmap, map<string, bool>* rcmap);
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << endl;
        cout << "   Function: given a genetic map file and fasta: "          << endl;
        cout << "             concatenate contigs into chrs. "               << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "   (version 1.0 - compiled on " << buildString << ")"                                   << endl << endl;
        cout << "   Usage: scaffolder upd_map_group1.txt contigs.fa "        << endl;
        cout << endl;
        exit(1);
    }
    double startT = clock();
    //
    string gmfile = (string)argv[1];
    vector<string> gmap;
    map<string, bool> rcmap; 
    if(!read_genetic_map(gmfile, &gmap, &rcmap))
    {
        return 1;
    }
    // fafile=manually_curated.fasta
    string fafile = (string)argv[2];
    ifstream ifp;
    ifp.open(fafile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << fafile << endl;
        return 1;
    }
    // traverse input fasta file
    map<string, string> seqs;
    int total_seq_num = 0;
    std::string line("");
    getline(ifp, line);    
    while(ifp.good())
    {       
        if(line.find(">") != std::string::npos)
        {
            string contigid = line.substr(1);
            cout << "   Info: checking " << total_seq_num+1 << "th seq: " << contigid << endl;                   
            
            if(std::find(gmap.begin(), gmap.end(), contigid) != gmap.end() ||
               std::find(gmap.begin(), gmap.end(), contigid+"sA") != gmap.end() ||
               std::find(gmap.begin(), gmap.end(), contigid+"sB") != gmap.end())
            {            
                string seq_name = contigid;
                total_seq_num ++;
                unsigned long chrlen = 0;
                string seq("");
                int seq_line_id = 0;    
                while(ifp.good())
                {
                    std::string chrseq("");                
                    getline(ifp, chrseq);
                    seq_line_id ++;
                    line    = chrseq;
                    if(line.find(">")!=std::string::npos) break; // next seq 
                    chrlen += chrseq.size();
                    if(seq_line_id > 1) seq += '\n';
                    seq    += chrseq;
                }
                seqs.insert(std::pair<string, string>(seq_name, seq));
            }else
            {
                getline(ifp, line);                
            }
        }else
        {
            getline(ifp, line);
        }
    }
    ifp.close();
    //
    vector<string> gmfileinfo = split_string(gmfile, '/');
    string gmflag  = gmfileinfo[gmfileinfo.size()-1];
    string outfile = "./scaffolded_" + gmflag.substr(0, gmflag.size()-4) + ".fa";
    ofstream ofp;
    ofp.open(outfile.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file " << outfile << endl;
        return 1;
    }
    ofp << ">" << gmflag.substr(0, gmflag.size()-4) << endl;
    //
    vector<string>::iterator ctgitr;
    vector<string>::iterator ctgitr_end;
    ctgitr     = gmap.begin();
    ctgitr_end = gmap.end();
    bool firstid = true;
    while(ctgitr != ctgitr_end)
    {
        string ctgid2 = *ctgitr;
        string ctgid  = *ctgitr;
        unsigned long separatepos = 0;
        // update special cases with contig id changed in marker list
        if(ctgid2.find("tig00003826_pilonsA")!=std::string::npos)
        {
            separatepos = 220000; // 210579
            ctgid = "tig00003826_pilon";
        }else
        if(ctgid2.find("tig00003826_pilonsB")!=std::string::npos)
        {
            separatepos = 220000; // 210579            
            ctgid = "tig00003826_pilon";            
        }else
        if(ctgid2.find("tig00003688_pilonsA")!=std::string::npos)
        {
            separatepos=2050000; 
            ctgid = "tig00003688_pilon";                              
        }else 
        if(ctgid2.find("tig00003688_pilonsB")!=std::string::npos)
        {
            separatepos=2050000; 
            ctgid = "tig00003688_pilon";                              
        }else ;        
        bool   rc    = rcmap[ctgid];
        map<string, string>::iterator seqitr = seqs.find(ctgid);
        if(seqitr == seqs.end())
        {
           cout << "   Error: inconsistent fasta file: " << ctgid << " not found. " << endl;
           return 1;
        }
        //
        if(!firstid)
        {
            ofp << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                << endl
                << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                << endl;            
        }      
        //ofp << ctgid << endl;  
        string this_seq2 = (*seqitr).second;
        string this_seq  = this_seq2;
        if(ctgid2.find("tig00003826_pilonsA")!=std::string::npos)
        {
            separatepos = 220000; // 210579
            this_seq = this_seq2.substr(0, separatepos);
        }else
        if(ctgid2.find("tig00003826_pilonsB")!=std::string::npos)
        {
            separatepos = 220000; // 210579            
            this_seq = this_seq2.substr(separatepos);
        }else
        if(ctgid2.find("tig00003688_pilonsA")!=std::string::npos)
        {
            separatepos=2050000; 
            this_seq = this_seq2.substr(0, separatepos);
        }else 
        if(ctgid2.find("tig00003688_pilonsB")!=std::string::npos)
        {
            separatepos=2050000; 
            this_seq = this_seq2.substr(separatepos);
        }else ;          
        if(rc)
        {
            string rcseq;
            reverse_complement(&this_seq, &rcseq);
            ofp << rcseq << endl;
            //ofp << this_seq << endl;                        
        }else
        {        
            ofp << this_seq << endl;
        }
        if(firstid)
        {
            firstid = false;
        }        
        //
        ctgitr ++;
    }
        
    ofp.close();
    //
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return 0;
}
//
bool read_genetic_map(string gmfile, vector<string>* gmap, map<string, bool>* rcmap)
{
    ifstream ifp;
    ifp.open(gmfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gmfile << endl;
        return false;
    }
    while(ifp.good())
    { 
        string line1("");
        getline(ifp, line1);
        if(line1.size()==0) continue;
        if(line1[0] == '#') continue;        
        vector<string> lineinfo1 = split_string(line1, '\t');       
         
        string line2("");
        if(!ifp.good())
        {
            cout << "   Error: not paired marker for " << line1 << endl;
            return false;
        }
        getline(ifp, line2); 
        if(line2.size()==0)
        {
            cout << "   Error: not paired marker for " << line1 << endl;
            return false;            
        }               
        vector<string> lineinfo2 = split_string(line2, '\t');        
        if(lineinfo1[0].compare(lineinfo2[0]) != 0)
        {
            cout << "   Error: not paired marker for " <<  line1 << endl;
            return false;
        }else
        {
            bool rc = (line1.find("right")!=std::string::npos);
            (*gmap).push_back(lineinfo1[0]);
            (*rcmap).insert(std::pair<string, bool>(lineinfo1[0], rc));
            cout << "   Info: " << lineinfo1[0] << " rev&comp = " << rc << endl;
        }
    }    
    ifp.close();
    cout << "   Info: " << (*gmap).size() << " contigs in the genetic map." << endl;
    return true;
}
//
bool reverse_complement(string* kmer, string* rckmer)
{
    if((*kmer).size()==0) return false;
    //
    *rckmer = *kmer;
    /* reverse sequence    */
    std::reverse((*rckmer).begin(), (*rckmer).end());
    /* lowercase sequence  */
    std::transform((*rckmer).begin(), (*rckmer).end(), (*rckmer).begin(), ::tolower);
    /* complement sequence */
    std::replace((*rckmer).begin(), (*rckmer).end(), 'a', 'T');
    std::replace((*rckmer).begin(), (*rckmer).end(), 'c', 'G');
    std::replace((*rckmer).begin(), (*rckmer).end(), 'g', 'C');
    std::replace((*rckmer).begin(), (*rckmer).end(), 't', 'A');
    //
    return true;
}
