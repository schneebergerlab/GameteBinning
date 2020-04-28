/*
    Given 
        a list of read counts (by "bedtools coverage -counts") at del_marker_finder-defined del-like markers
    1. normalize read counts according to total aligned reads        
    2. define/classify hap-chr and hom-chr regions.
    2019-12-09 started
    2020-01-08 updated: using RPKM for genotyping a del-like region.
    Written by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
*/
#include            <map>
#include         <string>
#include         <vector>
#include        <fstream>
#include        <sstream>
#include       <iostream>
#include      <algorithm>
#include       <string.h>
#include       <stdlib.h>
#include       <assert.h>
#include     <sys/stat.h>
#include       <dirent.h>
#include         <math.h>
#include "split_string.h"
using namespace std;
// del-like marker
struct DNODE
{
    string        ctg;
    unsigned long sta;
    unsigned long end;    // end of del-like marker
    unsigned long dep;    // average depth determined by leaf illumina sequencing     
    string        haphom; // class: hap_high;hap_low;hom_low;hom_high
};
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fpollen,
                 string* fdepth,
                 int*    colsample,
                 int*    colbarcode,
                 string* outprefix,
                 bool*   verbose); 
bool get_pollens(string  fpollen, 
                 int     colsample,
                 int     colbarcode,
                 map<string, string>* pollens,
                 map<int, string>*    pollens_ids,
                 map<string, int>*    pollens_ids_rev);
bool get_leaf_depth(string                                 fdepth, 
                 map<string, map<unsigned long, DNODE> >*  leafdepth);
bool get_read_stat(map<string, string>                     pollens,
                 map<string, int>                          pollens_ids_rev,
                 string                                    tmpfolder_s0,
                 map<string, map<unsigned long, DNODE> >   leafdepth,
                 map<string, map<string, unsigned long> >* all_pollen_read_counts,
                 map<string, map<string, unsigned long> >* all_pollen_read_counts_rpkm);
bool get_aligned_readnum(string this_align_file, 
                 unsigned long* align_raw, 
                 double*        align_rate, 
                 unsigned long* align_num);
bool get_region_readnum(string               this_depth_file, 
                 map<string, unsigned long>* this_pollen_read_counts);
bool genotype_del_like_marker(map<string, map<unsigned long, DNODE> >  leafdepth, 
                              string                                   tmpfolder_s, 
                              map<int, string>                         pollens_ids,
                              map<string, map<string, unsigned long> > all_pollen_read_counts,
                              map<string, map<string, unsigned long> > all_pollen_read_counts_rpkm);   
string get_depth_code(unsigned long this_depth);                                            
//
int main(int argc, char* argv[])
{
    if(argc < 11)
    {
        // g++ del_marker_genotyper.cpp split_string.cpp -O3 -o del_marker_genotyper
        cout << "\n   Function: genotype \"deletion-like\" marker (=regions without snp markers). ";
        const char *buildString = __DATE__ ", " __TIME__ ".";        
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "             Note: \"deletion\" can be: real or artifacial resulted from the case that " << endl
             << "                   a homologous haplotig diverged from the other thus reads from the other unmapped." << endl
             << "                   Markers defined with del_marker_finder; reads counts collected with bedtools coverage -counts."
             << endl;
        cout << "\n   Usage: "
             << "del_marker_genotyper --pollen pollen_consensus_flist.txt --leaf-depth del_like_interval_avg_depth_sorted.bed --sample 8 --barcode 9 -o prefix_out" << endl << endl 
             << "          Note: --sample 8 means sample info at 8th cell of splitted path, while barcode at 9th. "  << endl
             << "          Note: --leaf-depth is the resulted file from del_depth_finder."                           << endl
             << "          Note: read counts will be searched under folders of cells in pollen_consensus_flist.txt:" << endl
             << "                for example, target paired-files include: " << endl
             << "                    path/sample_A/AAACGGGTCCAGGCAT/AAACGGGTCCAGGCAT_del_like_read_count, and "      << endl
             << "                    path/sample_A/AAACGGGTCCAGGCAT/bowtie2.err"
             << endl
             << endl;
        return 1;
    }
    double startT= clock();
    // step 0. initialize from the cmd line
    cout << "\n   Initialize from cmd line: " << endl << endl;
    string fpollen   = "";
    string fdepth    = "";    
    int    colsample = 8;
    int    colbarcode= 9;
    string outprefix = "fun_";
    bool   verbose   = false; 
    if(!get_options(argc,
                    argv,
                    &fpollen,
                    &fdepth,
                    &colsample,
                    &colbarcode,
                    &outprefix,
                    &verbose) )
    {
        return 1;
    }
    // step 1. read pollen cell info
    map<string, string> pollens;        // <path, barcode>
    map<int, string>    pollens_ids;    // <pollenid, path> // 1-based id: for example, [1, 895] in apricot
    map<string, int>    pollens_ids_rev;// <path, pollenid> // 1-based id: 
    if(!get_pollens(fpollen, 
                    colsample,
                    colbarcode,
                    &pollens,
                    &pollens_ids,
                    &pollens_ids_rev)   )
    {
        cout << "   Error: cannot read pollen info. " << endl;
        return false;
    }
    // step 2. read leaf depth info at del-like regions.
    map<string, map<unsigned long, DNODE> > leafdepth;
    if(!get_leaf_depth(fdepth, &leafdepth) )
    {
        cout << "   Error: cannot read leaf-depth info at del-like regions. " << endl;
        return false;
    }                     
    //
    string tmpfolder_s = outprefix+"_tmp_pollen_del_like_genotypes_addi";
    DIR* dir = opendir(tmpfolder_s.c_str());
    if(dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder_s.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder_s << endl;
            return false;
        }
    }
    else;
    //    
    string tmpfolder_s0 = tmpfolder_s + "/s0_genotype_pollen_seq_del_wise";
    dir = opendir(tmpfolder_s0.c_str());
    if(dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder_s0.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder_s0 << endl;
            return false;
        }
    }
    else;    
    // step 3. for each pollen, get read count at all del-like regions
    //         result: string ofilename = tmpfolder_s0 + "/" + sample_A/B + "_" + barcode + "_counts_normalized.bed";
    map<string, map<string, unsigned long> > all_pollen_read_counts;
    map<string, map<string, unsigned long> > all_pollen_read_counts_rpkm;
    if(!get_read_stat(pollens, 
                      pollens_ids_rev, 
                      tmpfolder_s0, 
                      leafdepth, 
                      &all_pollen_read_counts, 
                      &all_pollen_read_counts_rpkm))
    {
        cout << "   Error: cannot read depths. " << endl;
        return false;
    }    
    // step 4. for each each del-like marker, integrate read counts of all pollens and derive "PM" pattern.
    if(!genotype_del_like_marker(leafdepth, 
                                 tmpfolder_s, 
                                 pollens_ids, 
                                 all_pollen_read_counts,
                                 all_pollen_read_counts_rpkm))
    {
        cout << "   Error: genotype del-like markers failed. " << endl;
        return false;
    }
    
    
    
    
    
    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool genotype_del_like_marker(map<string, map<unsigned long, DNODE> >  leafdepth, 
                              string                                   tmpfolder_s, 
                              map<int, string>                         pollens_ids,
                              map<string, map<string, unsigned long> > all_pollen_read_counts,
                              map<string, map<string, unsigned long> > all_pollen_read_counts_rpkm)
{
    /*
        INPUT:        
           leafdepth     gives the contig-position-sorted del-like regions. 
                         <ctg, <pos, {ctg,sta,end,dep,hap-hom} > >
           tmpfolder_s   gives the folder to collect resulted files
           pollens_ids   gives order of pollens from original 895 input cell list, and path/to/sample_x/barcode
                         <1-895, path/to/sample/barcode> 
                         this order of 1-895 needs to be preserved to match SNP-based PM pattern along a contig later.                           
           all_pollen_read_counts, gives raw read counts for all pollens (<= RPKM used yet -- they are effective)                     
                         <path/to/sample/barcode, <ctg\tsta\tend, raw_depth> >   
           all_pollen_read_counts_rpkm, gives normalized read counts ...
                         <path/to/sample/barcode, <ctg\tsta\tend, rpkm_depth> >  
                         Example: if matching the P pattern below to snp-marker contigs, 
                                  errors (given by *) tend to occur with raw_cnt definition, but not with rpkm definition.

--from del - tig00001216_pilon   65857   79333                         
206c024c2c22a0aad1000002c4000d40ac200ca0100ab000014b03000000b8000b0000a220a00002b8a800000000c20a2000000b04e000000a0000000000d00000a1000ae3b00a320a00030ec16e00010000000002a0020100 raw
a0b40a350401508540000006460007b0445004300003300004b7060000003b0007000044503000085a5a0000000043023000000409600000030000000000500000310002503002000600000847c70001000000000530090000 rpkm
PUPPUPPPPPPPPUPPPPUUUUUPPPUUUPPUPPPUUPPUPUUPPUUUUPPPUPUUUUUUPPUUUPUUUUPPPUPUUUUPPPPPUUUUUUUUPPUPPUUUUUUPUPPUUUUUUPUUUUUUUUUUPUUUUUPPUUUPPPPUUPPPUPUUUPUPPPPPUUUPUUUUUUUUUPPUUPUPUU PM-del
|.||..||*|*||.|||*.....|||...||.|||..||.*..||....|||.|......||...|....|||.|....|||||........||.||......|.||......|..........|.....|*...||*|..|**.|...*.|||||...*.........||..|.*..
PPPPPUPPMPMPPMPPPMMMMMUPPPPPPPPPPPPPPPPMMMMPPMMMPPPPPPPMMMMPPPMMMPPUPPPPPPPPMPPPPPPPUMMMMMMPPPPPPMMMMMPPPPPPMMMMUPMMMMMMMMMMPUPMMMPMMMMPPMPPUPMMMPPMMMMPPPPPMMMMMMMMMPMMMPPPPPMMPP PM-snp
-- from snp - tig00003827_pilon   right

                                                
        OUTPUT:
           1. PM patterns for hap-del-like regions
              -- this is necessary:
                    if the contig is     in linkage group, find out which genotype should have pacbio reads assigned
                    if the contig is not in linkage group, insert it to linkage group according to PM pattern, and 
                                                           find out which genotype should have pacbio reads assigned
           2. PM patterns for hom-del-like regions 
              -- this is not necessary, as if it is hom, pacbio reads go to both haplotypes.                 
    */
    std::stringstream gtinfo;
    gtinfo << tmpfolder_s
           << "/s2_genotype_contig_seq_del_like.txt\0"; //
    ofstream ctgPMofp;
    ctgPMofp.open((gtinfo.str()).c_str(), ios::out);
    if(!ctgPMofp.good())
    {
        cout   << "   Error: cannot open file " << gtinfo.str() << endl;
        return false;
    }    
    // traverse all del-like markers
    map<string, map<unsigned long, DNODE> >::iterator oditr;
    map<string, map<unsigned long, DNODE> >::iterator oditr_end;
    oditr     = leafdepth.begin();
    oditr_end = leafdepth.end();
    while(oditr != oditr_end)
    {
        map<unsigned long, DNODE> this_del_like = (*oditr).second;
        map<unsigned long, DNODE>::iterator delitr;
        map<unsigned long, DNODE>::iterator delitr_end;
        delitr     = this_del_like.begin();
        delitr_end = this_del_like.end();
        //
        while(delitr != delitr_end)
        {
            DNODE this_node   = (*delitr).second; // current del-like region == current marker
            string        ctg = this_node.ctg;
            unsigned long sta = this_node.sta;
            unsigned long end = this_node.end;
            unsigned long dep = this_node.dep;    // raw
            string     haphom = this_node.haphom;
            //
            std::stringstream this_key;
            this_key.str("");
            this_key << ctg << "\t" << sta << "\t" << end;
            // traverse all pollens
            map<int, string>::iterator piditr;
            map<int, string>::iterator piditr_end;
            piditr     = pollens_ids.begin();
            piditr_end = pollens_ids.end();
            std::stringstream read_cnt_seq;
            read_cnt_seq.str("");
            std::stringstream read_cnt_seq_rpkm;
            read_cnt_seq_rpkm.str("");            
            std::stringstream PM_seq;
            PM_seq.str("");
            int cnt1 = 0; // how many effective pollen
            while(piditr != piditr_end)
            {
                string ipath = (*piditr).second; // path/sample/barcode
                // get raw read counts of this pollen
                assert(all_pollen_read_counts.find(ipath) != all_pollen_read_counts.end());
                map<string, unsigned long> this_pollen_read_counts = all_pollen_read_counts[ipath];   
                assert(this_pollen_read_counts.find(this_key.str()) != this_pollen_read_counts.end());
                // get rpkm read counts of this pollen
                assert(all_pollen_read_counts_rpkm.find(ipath) != all_pollen_read_counts_rpkm.end());
                map<string, unsigned long> this_pollen_read_counts_rpkm = all_pollen_read_counts_rpkm[ipath];   
                assert(this_pollen_read_counts_rpkm.find(this_key.str()) != this_pollen_read_counts_rpkm.end());                
                // abcdefghij
                read_cnt_seq      << get_depth_code(this_pollen_read_counts[this_key.str()]);        
                read_cnt_seq_rpkm << get_depth_code(this_pollen_read_counts_rpkm[this_key.str()]);                        
                //
                if(this_pollen_read_counts_rpkm[this_key.str()] > 0) // exist
                {
                    PM_seq       << "P"; // note: for hom regions, "P"=="U" <= equal probability of being "P" and "M"
                    cnt1 ++;
                }else
                {
                    PM_seq       << "U";                  
                }
                piditr ++;
            }
            // output this del-like marker
            ctgPMofp << read_cnt_seq.str()      << "\t" << this_key.str() << "\t" << haphom << "\t" << cnt1 << "/895\t" 
                     << "raw_cnt" << endl;
            ctgPMofp << read_cnt_seq_rpkm.str() << "\t" << this_key.str() << "\t" << haphom << "\t" << cnt1 << "/895\t" 
                     << "rpkm"    << endl;            
            ctgPMofp << PM_seq.str()            << "\t" << this_key.str() << "\t" << haphom << "\t" << cnt1 << "/895\t" 
                     << "pmpat"   << endl;
            //
            delitr ++;
        }
        //
        ctgPMofp << endl;                              
        //
        oditr ++;
    }
    //
    ctgPMofp.close();
    //
    return true;
}
string get_depth_code(unsigned long this_depth)
{
    std::stringstream read_cnt_seq;
    read_cnt_seq.str("");
    if(this_depth < 10)
    {
        // [0,10)
        read_cnt_seq << this_depth;
    }else
    if(this_depth >= 10 &&
       this_depth <  20)           
    {
        read_cnt_seq << "a";
    }else
    if(this_depth >= 20 &&
       this_depth <  30)           
    {
        read_cnt_seq << "b";
    }else
    if(this_depth >= 30 &&
       this_depth <  40)           
    {
        read_cnt_seq << "c";
    }else
    if(this_depth >= 40 &&
       this_depth <  50)           
    {
        read_cnt_seq << "d";
    }else
    if(this_depth >= 50 &&
       this_depth <  60)           
    {
        read_cnt_seq << "e";
    }else
    if(this_depth >= 60 &&
       this_depth <  70)           
    {
        read_cnt_seq << "f";
    }else
    if(this_depth >= 70 &&
       this_depth <  80)           
    {
        read_cnt_seq << "g";
    }else
    if(this_depth >= 80 &&
       this_depth <  90)           
    {
        read_cnt_seq << "h";
    }else
    if(this_depth >= 90 &&
       this_depth <  100)           
    {
        read_cnt_seq << "i";
    }else 
    {
        // (100, inf)
        read_cnt_seq << "j";
    }   
    // 
    return read_cnt_seq.str();
}
//
bool get_read_stat(map<string, string> pollens, 
                   map<string, int>    pollens_ids_rev,
                   string              tmpfolder_s0, 
                   map<string, map<unsigned long, DNODE> >  leafdepth,
                   map<string, map<string, unsigned long> >* all_pollen_read_counts,
                   map<string, map<string, unsigned long> >* all_pollen_read_counts_rpkm)
{
    /*  
        RPKM = numReads / ( delLikeLength/1000 * totalNumReads/1000000)
        CPM  = numReads * 1.0 / totalNumReads * 1000000
    */
    if(pollens.size()==0)
    {
        cout << "   Error: no pollen data found. " << endl;
        return false;
    }
    //
    map<string, string>::iterator piditr;
    map<string, string>::iterator piditr_end;
    
    piditr     = pollens.begin();
    piditr_end = pollens.end();
    while(piditr != piditr_end)
    {
        string ipath = (*piditr).first; // path including barcode info
        string ibc   = (*piditr).second;// barcode 
        vector<string> ipathinfo = split_string(ipath, '/');        
        // total reads aligned
        string this_align_file  = ipath + "/" + "bowtie2.err";   
        unsigned long align_raw;
        double        align_rate; 
        unsigned long align_num = 1; // avoid 0          
        if(!get_aligned_readnum(this_align_file, 
                                &align_raw, 
                                &align_rate, 
                                &align_num) )
        {
            cout << "   Error: failed in reading alignment info file: " << this_align_file << endl;
            return false;
        }        
        // del-regional reads aligned
        string depthfile = ipath + "/" + ibc + "_del_like_read_count.bed";
        map<string, unsigned long> this_pollen_read_counts;
        map<string, unsigned long> this_pollen_read_counts_rpkm;
        if(!get_region_readnum(depthfile, &this_pollen_read_counts))
        {
            cout << "   Error: failed in reading depth file: " << depthfile << endl;
            return false;
        }
        this_pollen_read_counts_rpkm = this_pollen_read_counts; // to be normalized later
        // collect raw read counts
        if((*all_pollen_read_counts).find(ipath) != (*all_pollen_read_counts).end())
        {
           cout << "   Error: pollen read count data found twice: " << ipath << endl;
           return false;
        }else
        {
            (*all_pollen_read_counts).insert(std::pair<string, map<string, unsigned long> >(ipath, this_pollen_read_counts) );
        }
        // calculate RPKM
        string ofilename = tmpfolder_s0 + "/" + ipathinfo[ipathinfo.size()-2] + "_" + ibc + "_counts_normalized.bed";
        ofstream ofp;
        ofp.open(ofilename.c_str(), ios::out);
        if(!ofp.good())
        {
            cout << "   Error: cannot write output to file " << ofilename << endl;
            return false;
        }
        ofp << "#0Total_raw_readpair * Align_rate = Aligned_readpair  ::::  " 
            << align_raw  << " * "    
            << align_rate << " = "
            << align_num  << endl;
        ofp << "#0Contig\tStart\tEnd\tRaw_count\tNormalized_RPKM\tNormalized_CPM\t1_based_pollen_id=input_cell_order\tsample/barcode" << endl;
        //
        /*
        map<string, unsigned long>::iterator depth_itr;
        map<string, unsigned long>::iterator depth_itr_end;
        depth_itr     = this_pollen_read_counts.begin();
        depth_itr_end = this_pollen_read_counts.end();
        while(depth_itr != depth_itr_end)
        {
            string this_key        = (*depth_itr).first; // ctg \t sta \t end
            vector<string> keyinfo = split_string(this_key, '\t');
            unsigned long this_sta = strtoul(keyinfo[1].c_str(), NULL, 0);
            unsigned long this_end = strtoul(keyinfo[2].c_str(), NULL, 0);   
            unsigned long this_size= this_end - this_sta + 1;         
            //
            unsigned long normalized_RPKM = round( (*depth_itr).second * 1.0 / this_size * 1000 / (2*align_num) *1000000 );
            unsigned long normalized_CPM  = round( (*depth_itr).second * 1.0                    / (2*align_num) *1000000 );
            //
            ofp << this_key            << "\t" 
                << (*depth_itr).second << "\t"            
                << normalized_RPKM     << "\t" 
                << normalized_CPM      << endl;
            //
            depth_itr ++;
        }
        */       
        map<string, map<unsigned long, DNODE> >::iterator oditr;
        map<string, map<unsigned long, DNODE> >::iterator oditr_end;
        oditr     = leafdepth.begin();
        oditr_end = leafdepth.end();
        while(oditr != oditr_end)
        {
            map<unsigned long, DNODE> this_del_like = (*oditr).second;
            map<unsigned long, DNODE>::iterator delitr;
            map<unsigned long, DNODE>::iterator delitr_end;
            delitr     = this_del_like.begin();
            delitr_end = this_del_like.end();
            while(delitr != delitr_end)
            {
                DNODE this_node = (*delitr).second;                
                std::stringstream this_key;
                this_key.str("");
                this_key << this_node.ctg << "\t" << this_node.sta << "\t" << this_node.end;                
                map<string, unsigned long>::iterator depth_itr = this_pollen_read_counts.find(this_key.str());
                if(depth_itr == this_pollen_read_counts.end())
                {
                    cout << "   Error: " << this_key.str()   << " not found in pollen read counts. " << endl;
                    cout << "          current pollen file " << depthfile << endl;
                    cout << "   Error: del-like leaf-depth and pollen read count have different del-like regions. "
                         << endl;
                    return false;
                }
                unsigned long this_size= this_node.end - this_node.sta + 1;                         
                unsigned long normalized_RPKM = round( (*depth_itr).second * 1.0 / 
                                                       this_size * 1000          / 
                                                       (2*align_num) *1000000 );
                unsigned long normalized_CPM  = round( (*depth_itr).second * 1.0 / 
                                                       (2*align_num) *1000000 );                                                                         
                assert(pollens_ids_rev.find(ipath)!=pollens_ids_rev.end());
                ofp << this_key.str()         << "\t" 
                    << (*depth_itr).second    << "\t"            
                    << normalized_RPKM        << "\t" 
                    << normalized_CPM         << "\t"
                    << pollens_ids_rev[ipath] << "\t"
                    << ipathinfo[ipathinfo.size()-2] + "/" + ibc 
                    << endl;
                //
                this_pollen_read_counts_rpkm[this_key.str()] = normalized_RPKM;
                //                                    
                delitr ++;
            }
            //
            oditr ++;
        }
        // collect normalized read counts
        if((*all_pollen_read_counts_rpkm).find(ipath) != (*all_pollen_read_counts_rpkm).end())
        {
           cout << "   Error: pollen read count data found twice: " << ipath << endl;
           return false;
        }else
        {
            (*all_pollen_read_counts_rpkm).insert(std::pair<string, map<string, unsigned long> >(ipath, this_pollen_read_counts_rpkm) );
        }        
        //
        ofp.close();
        //
        piditr ++;
    }    
    return true;
}
//
bool get_region_readnum(string this_depth_file, map<string, unsigned long>* this_pollen_read_counts)
{
    ifstream ifp;
    ifp.open(this_depth_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot find file " << this_depth_file << endl;
        return false;
    }    
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() ==0) continue;
        if(line[0]  == '#') continue;
        
        /*
            tig00000026_pilon	1	13743	2
            tig00000026_pilon	13743	37715	0
            tig00000026_pilon	37715	105094	2
            tig00000026_pilon	105094	173329	6        
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<4)
        { 
            cout << "   Warning: skipping insufficient line info at " << line << endl;
            continue;
        } 
        string this_key = lineinfo[0] + "\t" + lineinfo[1] + "\t" + lineinfo[2];
        int this_count  = strtoul(lineinfo[3].c_str(), NULL, 0);        
        (*this_pollen_read_counts).insert(std::pair<string, unsigned long>(this_key, this_count));
        //
    }
    ifp.close();
    return true;
}
//
bool get_aligned_readnum(string         this_align_file, 
                         unsigned long* align_raw, 
                         double*        align_rate, 
                         unsigned long* align_num)
{
    *align_raw   = 0;
    *align_rate  = 100; // percentage
    *align_num   = 0;
    ifstream ifp;
    ifp.open(this_align_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot find file " << this_align_file << endl;
        cout << "          Hint: the file should have at least two lines: " << endl
             << "                xxxx reads; of these:"          << endl
             << "                ... "                           << endl
             << "                xx.xx\% overall alignment rate" 
             << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.find("reads; of these:") != std::string::npos)
        {
            vector<string> lineinfo = split_string(line, ' ');
            *align_raw              = strtoul(lineinfo[0].c_str(), NULL, 0);            
        }else
        if(line.find("overall alignment rate") != std::string::npos)
        {
            vector<string> lineinfo = split_string(line, ' ');
            *align_rate             = atof(lineinfo[0].substr(0, lineinfo[0].size()-1).c_str());            
        }
    }
    ifp.close();    
    //
    *align_num = round( (*align_raw) * (*align_rate) / 100);
    //
    return true;
}
//
bool get_leaf_depth(string fdepth, 
                 map<string, map<unsigned long, DNODE> >* leafdepth)
{
    ifstream ifp;
    ifp.open(fdepth.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: failed in reading " << fdepth << endl;
        return false;
    }
    int del_like_num = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        // tig00000013_pilon	5117	8635	219	hap_low
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5) 
        {
            cout << "   Warning: insufficient line info at " << line << endl;
            continue;
        }
        DNODE this_node;
        this_node.ctg    = lineinfo[0];
        this_node.sta    = strtoul(lineinfo[1].c_str(), NULL, 0);
        this_node.end    = strtoul(lineinfo[2].c_str(), NULL, 0);
        this_node.dep    = strtoul(lineinfo[3].c_str(), NULL, 0);
        this_node.haphom = lineinfo[4];
        //
        if((*leafdepth).find(this_node.ctg) == (*leafdepth).end())
        {
            map<unsigned long, DNODE> this_del_like;
            this_del_like.insert(std::pair<unsigned long, DNODE>(this_node.sta, this_node));
            (*leafdepth).insert(std::pair<string, map<unsigned long, DNODE> >(this_node.ctg, this_del_like));
        }else
        {
            (*leafdepth)[this_node.ctg].insert(std::pair<unsigned long, DNODE>(this_node.sta, this_node));
        }
        del_like_num ++;                
    }
    cout << "   Info: "        << del_like_num << " del-like regions collected from " 
         << (*leafdepth).size() << " contigs. " << endl;
    return true;
}                 
//
bool get_pollens(string  fpollen, 
                 int     colsample,
                 int     colbarcode,
                 map<string, string>* pollens,
                 map<int, string>* pollens_ids,
                 map<string, int>* pollens_ids_rev)
{
    /* extract info from lines below:
       /path/to/sample_A/AAACGGGTCCAGGCGA/shoremap_converted/extracted_consensus_0_update.txt
    */
    ifstream ifp;
    ifp.open(fpollen.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << fpollen << endl;
        return false;
    }
    int ipl = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        //
        vector<string> lineinfo = split_string(line, '/');
        if(lineinfo.size() < colbarcode)
        {
            cout << "   Warning: insufficient line info at " << line << endl;
            continue;
        }
        string path("");
        for(int ii=0; ii<colsample+1; ii++) // path = "/path/to/sample_A/AAACGGGTCCAGGCGA"
        {
            path += "/";
            path += lineinfo[ii];    
        }
        string barcode = lineinfo[colbarcode-1];
        (*pollens).insert(std::pair<string, string>(path, barcode));
        ipl ++;
        (*pollens_ids).insert(std::pair<int, string>(ipl, path));
        (*pollens_ids_rev).insert(std::pair<string, int>(path, ipl));  
        cout << "      Check: " << ipl << "th\t" << lineinfo[colsample-1] << "/" << barcode << endl;
    }
    ifp.close();
    cout << "   Info: " << (*pollens).size() << " pollen info collected. " << endl;
    return true;
}
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fpollen,
                 string* fdepth,
                 int*    colsample,
                 int*    colbarcode,                 
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
        if(optstr.compare("--pollen") == 0)
        {
            ic ++;
            *fpollen = (string)argv[ic];
            ifstream fp;
            fp.open((*fpollen).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: pollen cells file provided: "         
                     << *fpollen << endl;
            }else
            {
                cout << "   Error: cannot open pollen cells file "      
                     << *fpollen << endl;
                return false;
            }
        }else
        if(optstr.compare("--leaf-depth") == 0)
        {
            ic ++;
            *fdepth = (string)argv[ic];
            ifstream fp;
            fp.open((*fdepth).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: del-like depth file provided: "         
                     << *fdepth << endl;
            }else
            {
                cout << "   Error: cannot open del-like depth file "      
                     << *fdepth << endl;
                return false;
            }
        }        
        else
        if(optstr.compare("-o") == 0)
        {
            ic ++;
            *outprefix = (string)argv[ic];
            cout << "   Info: output files will be labeled with \""    
                 << *outprefix << "\"" << endl;
        }else
        if(optstr.compare("--sample") == 0)
        {
            ic ++;
            *colsample = atoi(argv[ic]);
            cout << "   Info: sample info at " << *colsample << "th splitted path cell. "
                 << endl;
        }else
        if(optstr.compare("--barcode") == 0)
        {
            ic ++;
            *colbarcode = atoi(argv[ic]);
            cout << "   Info: barcode info at " << *colbarcode << "th splitted path cell. "
                 << endl;
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

