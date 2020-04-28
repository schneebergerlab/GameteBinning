/*
    Given 
        a list of asPollinator-phased snp markers
        linkage-analysis correction on P/Ms,
        linkage-analysis phased del-like P/Ms,
        linkage group of snp-contig-PM contigs,
        a bam file with pacbio reads aligned,
    separate pacbio reads into 2*n clusters.
    Note: there is 500 bp cutoff to consider del-like region overlapping a pacbio read.
    2019-11-09 Started: 
    2020-01-08 Updated: consider del-like markers (already phased to JoinMap-phased snp-marker contigs by asCaffolder_v2).
    Written by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
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
#include    <unistd.h>
#include "split_string.h"
using namespace std;
//
struct PAMA
{
    string patallele;
    string matallele;
    string checkallele;
};
// del-like marker
struct DELMARKER
{
    string        ctg;
    unsigned long sta;
    unsigned long end;     // end of del-like marker
    string        haphom;  // class: hap_high;hap_low;hom_low;hom_high (if hom, reads go to PPP and MMM; otherwise go 'geno')
    string        geno;    // phased PM value/pmpat for this del-like marker; M indicates flipped already according to phased snp-contig-PMs
    string        gphase;  // belong to which phasing group
    string        gmap;    // belong to which genetic map
};
double minCOscore  = 0.64; // 
//
bool get_options(int                                     argc,
                 char*                                   argv[],
                 string*                                 fsam,
                 string*                                 fmarker,
                 string*                                 fmarkerdel,
                 string*                                 fphase,
                 string*                                 outprefix, 
                 bool*                                   verbose,
                 double*                                 minCOscore);
int decipher_cigar(string                                cigar, 
                 vector<char>*                           operation, 
                 vector<int>*                            count);
bool get_contig_phase_files(string                       fphase, 
                 vector<string>*                         groupphase);
bool get_contig_phase(string                             jmfile, 
                 map<string, bool>*                      ctgJoinMapPhasing,
                 map<string, bool>*                      ctgJoinMapPhasing_all,
                 map<string, string>*                    ctgGroup);  
bool get_del_like_phase(string                           fmarkerdel, 
                 map<string, map<unsigned long, DELMARKER> >* del_like_marker,
                 map<string, string>*                    ctgGroup);
bool read_asPollinator_phased_marker(string              fmarker,
                 map<string, bool>                       ctgJoinMapPhasing_all, 
                 map<string, map<unsigned long, PAMA> >* correctedmkr);
bool align_read_to_ref(vector<char>                      operation, 
                 vector<int>                             count, 
                 string*                                 orignal_read, 
                 string*                                 updated_read); 
bool genotype_pbread(string*                             updated_read, 
                 map<unsigned long, PAMA>*               this_marker, 
                 unsigned long                           start, 
                 unsigned long                           end, 
                 string*                                 geno,
                 map<unsigned long, PAMA>*               intersectMkr);
string get_break_pos(string                              poPat, 
                 string*                                 pmflag,
                 unsigned long*                          leftp, 
                 double*                                 score);   
bool overlap_read_del_like(map<unsigned long, DELMARKER>  this_del_like, 
                           unsigned long                  read_sta, 
                           unsigned long                  read_end,
                           map<unsigned long, DELMARKER>* this_del_like_overlapped,
                           string*                        del_phase);
//
int main(int argc, char* argv[])
{
    if(argc < 13)
    {
        // g++ pacbio_genotyper.cpp split_string.cpp -O3 -o pacbio_genotyper
        cout << "\n   Function: separate pacbio reads using genotypes according to phased contig-snp/del markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "pacbio_genotyper"
             << " --sam read_alignment.sam"
             << " --marker snp_list.txt"
             << " --marker2 del_like_list.txt"
             << " --phase ctg_phasing.txt"
             << " --ims 0.64 "
             << " -o prefix_out" << endl << endl;
        cout << "   Note: "      << endl;
        cout << "      --phase is a file listing n files of snp-contig-PM phased contigs, n=linkgage group number. " 
             << endl;
        cout << "      --ims   provides the minimum score to report a co along a read [0.51, 1]."
             << endl;
        cout << "      --marker2 is from asCaffolder_v2: flag_s2_genotype_contig_seq_del_like.txt. " << endl
             << endl;
        return 1;
    }
    double startT= clock();
    // step 0. initialize from the cmd line
    cout << "Initialize input from cmd line: " << endl << endl;
    string fsam      = "";
    string fmarker   = "";
    string fmarkerdel= "";    
    string fphase    = "";
    string outprefix = "fun_";
    bool   verbose   = false; 
    if(!get_options(argc,
                    argv,
                    &fsam,
                    &fmarker,
                    &fmarkerdel,
                    &fphase,
                    &outprefix, 
                    &verbose,
                    &minCOscore) )
    {
        return 1;
    } 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// step 1. get contig phasing: snp-contig markers and del-like markers.
    /// step 1.1 snp-marker based phasing/grouping info
    cout << endl;
    vector<string> groupphase;      
    if(!get_contig_phase_files(fphase, &groupphase))
    {
        return 1;
    }
    //
    map<string, map<string, bool> > contig_phasing;        // <this_linkage_file, <contig_id, bool_flip > >
    map<string, bool>               ctgJoinMapPhasing_all; // <contig_id, bool-flip>
    map<string, string>             ctgGroup;              // <contig_id, this_linksage_file>
    //
    cout << endl;
    vector<string>::iterator gpitr;
    vector<string>::iterator gpitr_end;
    gpitr     = groupphase.begin();
    gpitr_end = groupphase.end();
    while(gpitr != gpitr_end)
    {
        string this_jmfile = *gpitr; // this is phasing/grouping file - not genetic map file.
        //
        vector<string> fileinfo  = split_string(this_jmfile, '/');
        string this_linkage_file = fileinfo[fileinfo.size()-1];        
        cout << "   Info: read linkage-phasing info from " << fileinfo[fileinfo.size()-1] << endl;    
        //
        map<string, bool> ctgJoinMapPhasing;
        if(!get_contig_phase(this_jmfile, 
                             &ctgJoinMapPhasing,
                             &ctgJoinMapPhasing_all,
                             &ctgGroup))
        {
            cout << "   Error: failed in reading " << this_jmfile << endl;
            return 1;
        }
        //
        contig_phasing.insert(std::pair<string, map<string, bool> >(this_linkage_file, ctgJoinMapPhasing));
        cout << "         " << ctgJoinMapPhasing.size() << " contigs phasing info collected. " << endl;
        //
        gpitr ++;
    }
    cout << "   Info: overall, " << ctgJoinMapPhasing_all.size() << " contigs phasing info collected.\n" << endl; 
    // step 1.2 del-like marker based grouping info.   
    // update ctgGroup on phasing/grouping info with del-like markers; insert missing contigs into this variable.
    map<string, map<unsigned long, DELMARKER> > del_like_marker;//<contig_id, <start, {ctg,sta,end,haphom,geno,gphase,gmap}> >
    if(!get_del_like_phase(fmarkerdel, &del_like_marker, &ctgGroup))
    {
        cout << "   Error: cannot read del-like marker file. " << endl;
        return 1;
    }    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////// step 2. read markers with correction on PM alleles
    map<string, map<unsigned long, PAMA> > correctedmkr; // map<chr, map<position, {mat-allele, pat-allele} > >
    if(!read_asPollinator_phased_marker(fmarker,
                                        ctgJoinMapPhasing_all, 
                                        &correctedmkr) )
    {
        return false;
    }
    if(true)
    {
        vector<string> mkrfileinfo = split_string(fmarker, '/');
        string corrmarkerfile = "lgcorrected_" + mkrfileinfo[mkrfileinfo.size()-1];
        ofstream mkrofp;
        mkrofp.open(corrmarkerfile.c_str(), ios::out);
        if(!mkrofp.good())
        {
            cout << "   Error: cannot open file to write corrected markers: " << corrmarkerfile << endl;
            return 1;
        }else
        {
            cout << "   Info: correction on marker phasing with groups collected in " << corrmarkerfile << endl;
        }
        map<string, map<unsigned long, PAMA> >::iterator citr;
        map<string, map<unsigned long, PAMA> >::iterator citr_end;
        citr     = correctedmkr.begin();
        citr_end = correctedmkr.end();
        while(citr != citr_end)
        {
            string thisctg = (*citr).first;
            //
            map<string, bool>::iterator phaseitr = ctgJoinMapPhasing_all.find(thisctg);
            int phasevalue = -1; // nonchanged due to no phasing info
            if(phaseitr != ctgJoinMapPhasing_all.end())
            {
                phasevalue = (int)ctgJoinMapPhasing_all[thisctg]; // 1-true-swapped; 0-false-kept
            }
            //
            map<unsigned long, PAMA>::iterator positr;
            map<unsigned long, PAMA>::iterator positr_end;
            positr     = correctedmkr[thisctg].begin();
            positr_end = correctedmkr[thisctg].end();
            while(positr != positr_end)
            {
                unsigned long thispos = (*positr).first;
                mkrofp << thisctg << "\t"
                       << thispos << "\t"
                       << correctedmkr[thisctg][thispos].matallele << "\t"
                       << correctedmkr[thisctg][thispos].patallele << "\t"
                       << "flip=" << phasevalue                    << endl;
                positr ++;
            }
            citr ++;
        }
        mkrofp.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// step 3. genotype pacbio read in sam
    // output files
    // create an intermediate folder for collecting details about a CO-molecule
    string tmpfolder = outprefix + "_snp_marker_separated_pbreads";
    DIR* dir = opendir(tmpfolder.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder << endl;
            return false;
        }
    }
    else;  
    // file handles on output files: <infilename+"_PPP/MMM_pbreads.fa", *ofp >    
    map<string, ofstream*> allofp;         
    // collecting group-specific unique read names: <infilename+"_PPP/MMM_pbreads.fa", <pbreadname, cnt> >               
    map<string, map<string, int> > groupedpbreadnames;
    // no real data - just for initializing groupedpbreadnames.    
    map<string, int> init_tmp;     
    map<string, map<string, bool> >::iterator grpfileitr;
    map<string, map<string, bool> >::iterator grpfileitr_end;
    grpfileitr     = contig_phasing.begin();
    grpfileitr_end = contig_phasing.end();
    while(grpfileitr != grpfileitr_end)
    {
        string this_infilename = (*grpfileitr).first;
        //
        string this_outfilename = tmpfolder + "/" + this_infilename + "_MMM_pbreads.fa";
        ofstream *mmmofp = new ofstream(this_outfilename.c_str(), ios::out);
        if(!(*mmmofp).good())
        {
            cout << "   Error: cannot set up output file for " << this_outfilename << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ofstream*>(this_infilename + "_MMM_pbreads.fa", mmmofp));
        groupedpbreadnames.insert(std::pair<string, map<string, int> >(this_infilename + "_MMM_pbreads.fa", init_tmp));                
        //
        this_outfilename = tmpfolder + "/" + this_infilename + "_PPP_pbreads.fa";
        ofstream *pppofp = new ofstream(this_outfilename.c_str(), ios::out);
        if(!(*pppofp).good())
        {
            cout << "   Error: cannot set up output file for " << this_outfilename << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ofstream*>(this_infilename + "_PPP_pbreads.fa", pppofp));
        groupedpbreadnames.insert(std::pair<string, map<string, int> >(this_infilename + "_PPP_pbreads.fa", init_tmp));
        //
        grpfileitr ++;
    }
    // pacbio reads collected    
    map<string, int> pbreadnames;                    // all collected pb read names (unique of 2587118)
    map<string, vector<string> > pbreadnames_detail; // all collected pb read names (unique of 2587118) <readname, <grp, align_len> >   
    map<string, int> pbreadnames_nolg;     // all collected pb read names without linkage group info (unique of )
    map<string, int> pbreadnames_lg;       // all collected pb read names with    linkage group info and with    marker (unique of )
    map<string, int> pbreadnames_lg_nomkr; // all collected pb read names with    linkage group info and without marker (unique of )    
    // unmapped pacbio reads with cigar as "*"
    string unmappedfile   = tmpfolder + "/unmapped_starcigar_pb_reads.fa"; 
    ofstream *ofp         = new ofstream(unmappedfile.c_str(), ios::out);
    allofp.insert(std::pair<string, ofstream*>("unmapped_starcigar_pb_reads.fa", ofp));
    groupedpbreadnames.insert(std::pair<string, map<string, int> >("unmapped_starcigar_pb_reads.fa", init_tmp));    
    // mapped pacbio reads but no linkage group info
    string mappednolgfile = tmpfolder + "/mapped_no_lg_pb_reads.fa"; 
    ofstream *nolgofp     = new ofstream(mappednolgfile.c_str(), ios::out);
    allofp.insert(std::pair<string, ofstream*>("mapped_no_lg_pb_reads.fa", nolgofp));  
    groupedpbreadnames.insert(std::pair<string, map<string, int> >("mapped_no_lg_pb_reads.fa", init_tmp));          
    // pacbio reads showing potential crossovers
    string cofile         = tmpfolder + "/crossover_pb_reads.fa";     
    ofstream *coofp       = new ofstream(cofile.c_str(), ios::out);   
    allofp.insert(std::pair<string, ofstream*>("crossover_pb_reads.fa", coofp));
    groupedpbreadnames.insert(std::pair<string, map<string, int> >("crossover_pb_reads.fa", init_tmp));     
    // caution on unphased.txt 2020-01-08; we have one case in apricot, one contig tig00005012_pilon with no reads mapped.
             
    //
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";
    unsigned long numNoMA  = 0; // number of reads showing multiple alignment: only one kept.
    string exampleNoMA     = "";
    unsigned long numraw   = 0;    
    string exampleNoSEQ    = "";
    unsigned long numNoSEQ = 0;
    unsigned long lengrped = 0; // pb alignments clearly related to linkage groups - length 
    unsigned long grpedNo  = 0; // pb alignments clearly related to linkage groups; note 1 read may have >1 alignments
    unsigned long lennonlg = 0; // pb alignments not related to linkage groups - length       
    unsigned long nonlgNo  = 0; // pb alignments not related to linkage groups - number 
    unsigned long lennomkr = 0; // pb alignments related to linkage groups, however no marker exist can determine P/M cluster. - length
    unsigned long lgnomkr  = 0; // pb alignments related to linkage groups, however no marker exist can determine P/M cluster. - number 
    //
    ifstream samifp;
    samifp.open(fsam.c_str(), ios::in);
    if(!samifp.good())
    {
        cout << "   Error: cannot open sam file " << fsam << endl;
        return 1;
    }
    int srandsleep = 1;
    sleep(srandsleep);    
    while(samifp.good())
    {
        string line("");
        getline(samifp, line);
        if(line.size()==0) continue;
        if(line[0] == '@') continue;
        //
        numraw ++;
        if(numraw%100000 == 0)
        {
            cout << "   info: " << numraw << "th aligm..." << endl;
        }        
        //        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<11)
        {
            cout << "   Warning: insufficient line info, skipped: " << line << endl;
            continue;
        }
        string pbname    = lineinfo[0]; // read name
        string thisflag  = lineinfo[1]; // flag value
        string thisctg   = lineinfo[2]; // reference contig name; there are special case
        string thispos   = lineinfo[3]; // first matching reference coordinate
        string thiscigar = lineinfo[5]; // cigar string
        string thisseq   = lineinfo[9]; // read sequence         
        if(thisseq.compare("*") == 0) 
        {
            numNoSEQ ++;
            if(exampleNoSEQ.size()==0)
            {
                exampleNoSEQ = line;
                cout << endl
                     << "   Warning: there are alignments without clear sequence, secondary skipped, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue; // caution
        }
        // special case 1: cigar string as star: not aligned
        if(thiscigar.compare("*") == 0) 
        {
            // READ SET I
            (*allofp["unmapped_starcigar_pb_reads.fa"])  << ">" << pbname  << endl;
            (*allofp["unmapped_starcigar_pb_reads.fa"])  <<        thisseq << endl;
            //
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = line;
                cout << endl
                     << "   Warning: there are alignments without explicit CIGAR string, skipped., e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }            
            numNoCG ++;
            lenNoCG  += thisseq.size();
            //
            if(groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].find(pbname) == 
               groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].end())
            {
                groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames["unmapped_starcigar_pb_reads.fa"][pbname] += 1;
            }
            //     
            continue;
        }
        // special case 2:
        int hexflag = strtol(thisflag.c_str(), NULL, 0);
        if(hexflag > 255)  // do not skip such lines???????? to discuss!!!
        {
            numNoMA ++;
            if(exampleNoMA.size()==0)
            {
                exampleNoMA = line;
                cout << endl
                     << "   Warning: there are alignments being secondary/supplementary etc, skipped, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue;
        }
        // collect all informative reads
        if(pbreadnames.find(pbname) == pbreadnames.end())
        {
            pbreadnames.insert(std::pair<string, int>(pbname, 1));
        }else
        {
            pbreadnames[pbname] += 1;
        }
        cout << endl;
        string alignedtogrp = "";
        cout << "   Check: read: " << pbname << " aligned to " << thisctg << "\t" << thispos << endl;               
        // If it contains 0x10, SEQ in SAM has been reversed+complemented from read in fastq
        string reversed = "nrc";
        if((hexflag & 0x10) == 0x10) reversed = "rc";
        //
        vector<char>  operation;
        vector<int>   count;        
        // covered len by alignment
        int           covlen  = decipher_cigar(thiscigar, &operation, &count);   
        // alignment start and end position on ref seq: these would intersect interested markers for current pb read
        unsigned long firstRefsMatch = strtoul(thispos.c_str(), NULL, 0);                       
        unsigned long spansecond     = firstRefsMatch+covlen-1; // 1-based         
        cout << "          span: " << firstRefsMatch              << "-"   << spansecond 
             << " => "             << spansecond-firstRefsMatch+1 << " bp" << endl; // 1-based 
        // update special cases with contig id changed in marker list
        if(thisctg.find("tig00003826_pilon")!=std::string::npos)
        {
            unsigned long separatepos = 220000; // 210579
            if(spansecond < separatepos)
            {
                thisctg = "tig00003826_pilonsA";                 
            }else
            if(firstRefsMatch > separatepos)
            {
                thisctg = "tig00003826_pilonsB";
            }else
            if(separatepos-firstRefsMatch >= spansecond-separatepos)
            {
                thisctg = "tig00003826_pilonsA";
            }else
            {
                thisctg = "tig00003826_pilonsB";
            }            
            cout << "         Warning: special contig " << thisctg << endl;
        }else
        if(thisctg.find("tig00003688_pilon")!=std::string::npos)
        {
            unsigned long separatepos=2050000; 
            if(spansecond < separatepos)
            {
                thisctg = "tig00003688_pilonsA";
            }else
            if(firstRefsMatch > separatepos)
            {
                thisctg = "tig00003688_pilonsB";
            }else
            if(separatepos-firstRefsMatch >= spansecond-separatepos)
            {
                thisctg = "tig00003688_pilonsA";
            }else
            {
                thisctg = "tig00003688_pilonsB";
            }
            cout << "         Warning: special contig " << thisctg << endl;                
        }else ;        
        //
        string updated_read("");
        if(!align_read_to_ref(operation, count, &thisseq, &updated_read))
        {
            cout << "   Warning: unexpected read " << thisseq << endl;
            continue;
        }
        cout << "   Check: original read in sam: " 
             << thisseq.size()      << " bp; " << endl;
        cout << "   Check: updated read with \'I\' removed, and \'D:-\' added: " 
             << updated_read.size() << " bp. " << endl;
        // prepare all markers along this contig        
        map<unsigned long, PAMA> this_contig_marker;
        this_contig_marker.clear();
        map<string, map<unsigned long, PAMA> >::iterator cpitr = correctedmkr.find(thisctg);     
        if(cpitr != correctedmkr.end())
        {
            this_contig_marker = correctedmkr[thisctg];
        }
        cout << "   Check: there are " << this_contig_marker.size() << " snp markers for " << thisctg << endl;        
        // cout << updated_read << endl;
        // snp markers: get genotype of the read along intersected snp markers
        string                   thisgeno;
        map<unsigned long, PAMA> intersectMkr;
        genotype_pbread(&updated_read, 
                        &this_contig_marker, 
                        firstRefsMatch, 
                        spansecond, 
                        &thisgeno,
                        &intersectMkr);
        // del-like markers: for later, prepare variables related to del-like markers -- might be checked later
        map<unsigned long, DELMARKER> this_del_like_overlapped;    // not used yet
        this_del_like_overlapped.clear();
        string                        del_phase;                   // P M  or U
        // check crossover
        string pmflag(""); // P M or U
        unsigned long leftp;
        double        score;
        string pmstring = get_break_pos(thisgeno, &pmflag, &leftp, &score);
        // find the phase file
        map<string, string>::iterator c2jmitr = ctgGroup.find(thisctg); // ctgGroup <contig_id, this_linkage_file>
        string tmp_jmfile = "";
        if(c2jmitr != ctgGroup.end())
        {
            cout << "   Check: this read-contig is in linkage group: " << (*c2jmitr).second << endl;
            // note here including those pb reads with CO defined
            grpedNo ++;
            lengrped += thisseq.size();
            // pb reads in linkage group; in file name flag
            tmp_jmfile = (*c2jmitr).second;
            //
            if(pmflag.compare("P")==0)
            {
                // READ SET II: snp marker assigned with P genotype
                string tmp_jmfilep = tmp_jmfile + "_PPP_pbreads.fa";
                if(allofp.find(tmp_jmfilep) == allofp.end())
                {
                    cout << "   Warning skipping phasing group file " << tmp_jmfilep << " as it's not found. " << endl;
                    continue;
                }
                (*allofp[tmp_jmfilep]) << ">"   
                                       << pbname 
                                       << " +lg"
                                       << " span "         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond 
                                       << " snp-phased"                                       
                                       << endl;   // read name
                (*allofp[tmp_jmfilep]) << thisseq         
                                       << endl;   // read seq
                //
                if(groupedpbreadnames[tmp_jmfilep].find(pbname) == 
                   groupedpbreadnames[tmp_jmfilep].end())
                {
                    groupedpbreadnames[tmp_jmfilep].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[tmp_jmfilep][pbname] += 1;
                } 
                //
                alignedtogrp = tmp_jmfilep;               
            }else
            if(pmflag.compare("M")==0)
            {
                // READ SET III: snp marker assigned with M genotype
                string tmp_jmfilem = tmp_jmfile + "_MMM_pbreads.fa";  
                assert(allofp.find(tmp_jmfilem) != allofp.end());
                (*allofp[tmp_jmfilem])  << ">"   
                                       << pbname 
                                       << " +lg"
                                       << " span "         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond 
                                       << " snp-phased"                                        
                                       << endl;   // read name
                (*allofp[tmp_jmfilem]) << thisseq         
                                       << endl;   // read seq  
                //
                if(groupedpbreadnames[tmp_jmfilem].find(pbname) == 
                   groupedpbreadnames[tmp_jmfilem].end())
                {
                    groupedpbreadnames[tmp_jmfilem].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[tmp_jmfilem][pbname] += 1;
                } 
                //
                alignedtogrp = tmp_jmfilem;                                                                                 
            }else
            if(pmflag.compare("U")==0) // no phased snp marker overlapping the read
            {
                cout << "   Check: further checking on del-like markers: " << endl;
                //////////////////////////////////////////////////////////////////////////////////////////////////////
                // 20200108: in this case, we can check more on if the read overlaps del-like markers
                map<unsigned long, DELMARKER> this_del_like; 
                this_del_like.clear();
                map<string, map<unsigned long, DELMARKER> >::iterator dmitr = del_like_marker.find(thisctg);   
                if(dmitr != del_like_marker.end())
                {
                    this_del_like = del_like_marker[thisctg];                    
                    cout << "   Check: there are "   << this_del_like.size() 
                         << " del-like markers for " << thisctg 
                         << endl;
                    // find overlaping dels with current read span and get region-phase
                    overlap_read_del_like(this_del_like, 
                                          firstRefsMatch, 
                                          spansecond, 
                                          &this_del_like_overlapped,
                                          &del_phase);
                    if(del_phase.find("PP")!=std::string::npos) del_phase = "P"; else
                    if(del_phase.find("MM")!=std::string::npos) del_phase = "M";
                }else
                {
                    cout << "   Check: there are "   << this_del_like.size() 
                         << " del-like markers for " << thisctg 
                         << endl;                     
                    cout << "   Warning: unexpectedly no snp and no del markers found for " << thisctg << endl;
                    del_phase = "U";
                }
                // READ SET DEL_I: del-like marker assigned with P genotype
                srand (time(NULL));                                  
                double randra = rand()%1000000000*1.0 / 1000000000.0;                  
                if(del_phase.compare("P") == 0 || (del_phase.compare("U") == 0 && randra<=0.5))
                {
                    string tmp_jmfilep = tmp_jmfile + "_PPP_pbreads.fa";
                    if(allofp.find(tmp_jmfilep) == allofp.end())
                    {
                        cout << "   Warning skipping phasing group file " << tmp_jmfilep << " as it's not found. " << endl;
                        continue;
                    }
                    (*allofp[tmp_jmfilep]) << ">"   
                                           << pbname 
                                           << " +lg"
                                           << " span "         
                                           << thisctg 
                                           << ":"   
                                           << firstRefsMatch   
                                           << "-"          
                                           << spansecond
                                           << " del-like-phased" 
                                           << endl;   // read name
                    (*allofp[tmp_jmfilep]) << thisseq         
                                           << endl;   // read seq
                    //
                    if(groupedpbreadnames[tmp_jmfilep].find(pbname) == 
                       groupedpbreadnames[tmp_jmfilep].end())
                    {
                        groupedpbreadnames[tmp_jmfilep].insert(std::pair<string, int>(pbname, 1));
                    }else
                    {
                        groupedpbreadnames[tmp_jmfilep][pbname] += 1;
                    } 
                    //
                    alignedtogrp = tmp_jmfilep;                                             
                }
                // READ SET DEL_II: del-like marker assigned with M genotype
                if(del_phase.compare("M") == 0 || (del_phase.compare("U") == 0 && randra>0.5))
                {
                    string tmp_jmfilem = tmp_jmfile + "_MMM_pbreads.fa";
                    if(allofp.find(tmp_jmfilem) == allofp.end())
                    {
                        cout << "   Warning skipping phasing group file " << tmp_jmfilem << " as it's not found. " << endl;
                        continue;
                    }
                    (*allofp[tmp_jmfilem]) << ">"   
                                           << pbname 
                                           << " +lg"
                                           << " span "         
                                           << thisctg 
                                           << ":"   
                                           << firstRefsMatch   
                                           << "-"          
                                           << spansecond
                                           << " del-like-phased" 
                                           << endl;   // read name
                    (*allofp[tmp_jmfilem]) << thisseq         
                                           << endl;   // read seq
                    //
                    if(groupedpbreadnames[tmp_jmfilem].find(pbname) == 
                       groupedpbreadnames[tmp_jmfilem].end())
                    {
                        groupedpbreadnames[tmp_jmfilem].insert(std::pair<string, int>(pbname, 1));
                    }else
                    {
                        groupedpbreadnames[tmp_jmfilem][pbname] += 1;
                    }
                    //
                    alignedtogrp = tmp_jmfilem;                     
                }
                //
                //cout << "  check: del_phase=" << del_phase << " with randra = " <<  randra << endl;
                //if(alignedtogrp.size() == 0) cout << "   Warning: this grp name is null in del-check " << pbname << endl;
                //
                if(del_phase.compare("P")==0 || del_phase.compare("M")==0)
                {
                    // pb alignments related to linkage groups, del-like marker existing can determine P/M cluster.     
                    // with lg with mkr
                    if(pbreadnames_lg.find(pbname) == pbreadnames_lg.end())
                    {
                        pbreadnames_lg.insert(std::pair<string,int>(pbname, 1));
                    }
                }                 
                // under condition: pmflag.compare("U")==0
                if(del_phase.compare("U") == 0)
                {
                    // with lg without snp/del mkr
                    if(pbreadnames_lg_nomkr.find(pbname) == pbreadnames_lg_nomkr.end())
                    {
                        pbreadnames_lg_nomkr.insert(std::pair<string, int>(pbname, thisseq.size()));
                    }
                    // pb alignments is related to linkage groups, 
                    // however no del/snp marker exist can determine P/M cluster.     
                    // reads have been collected in both genotypes (see above)
                    lgnomkr  ++; 
                    lennomkr += thisseq.size();                           
                }
            }
            else ;            
            //
            if(pmflag.compare("P")==0 || pmflag.compare("M")==0)
            {
                // pb alignments related to phasing/linkage groups, SNP marker existing can determine P/M cluster.     
                // with lg with mkr
                if(pbreadnames_lg.find(pbname) == pbreadnames_lg.end())
                {
                    pbreadnames_lg.insert(std::pair<string,int>(pbname, 1));
                }
            } 
        }else
        {
            cout << "   Check: this read-contig is not in any linkage group " << endl;            
            cout << "          read would go to mapped_no_lg_pb_reads.fa. "   << endl; 
            nonlgNo  ++;
            lennonlg += thisseq.size();
            // without lg
            if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
            {
                pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
            }            
            // pb reads mapped but no linkage group info
            // READ SET IV            
            (*allofp["mapped_no_lg_pb_reads.fa"])  << ">"   << pbname << " no-lg "
                                                   << " span "        << thisctg      << ":"   
                                                   << firstRefsMatch  << "-"          << spansecond << endl;
            (*allofp["mapped_no_lg_pb_reads.fa"])  << thisseq         << endl; 
            //
            string tmp_jmfilenlg = "mapped_no_lg_pb_reads.fa";
            if(groupedpbreadnames[tmp_jmfilenlg].find(pbname) == 
               groupedpbreadnames[tmp_jmfilenlg].end())
            {
                groupedpbreadnames[tmp_jmfilenlg].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames[tmp_jmfilenlg][pbname] += 1;
            }
            //
            alignedtogrp = tmp_jmfilenlg;                                               
        }
        //
        if(leftp != 0)
        {
            // pb reads showing potential crossovers (note the reads can overlap READ SET II-III-IV)
            (*allofp["crossover_pb_reads.fa"])  << ">"   << pbname << " co-score " << score
                                                << " after "       << leftp        << "th/" << intersectMkr.size()
                                                << " mkr"
                                                << " span "        << thisctg      << ":"   
                                                << firstRefsMatch  << "-"          << spansecond << endl;
            (*allofp["crossover_pb_reads.fa"])  << thisseq         << endl;
            //
            string tmp_jmfileco = "crossover_pb_reads.fa";
            if(groupedpbreadnames[tmp_jmfileco].find(pbname) == 
               groupedpbreadnames[tmp_jmfileco].end())
            {
                groupedpbreadnames[tmp_jmfileco].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames[tmp_jmfileco][pbname] += 1;
            }                           
        }                            
        if(true && thisgeno.size()!=0)
        {
            cout << "   Check: pos mat-allele pat-allele this-allele genotyped ith-mkr-in-span" << endl;           
            map<unsigned long, PAMA>::iterator imitr;
            map<unsigned long, PAMA>::iterator imitr_end;
            imitr     = intersectMkr.begin();
            imitr_end = intersectMkr.end();
            int pi    = 0;
            while(imitr != imitr_end)
            {
                cout << "          " 
                     << (*imitr).first              << "\t"
                     << (*imitr).second.matallele   << "\t"
                     << (*imitr).second.patallele   << "\t"
                     << (*imitr).second.checkallele << "\t" 
                     << thisgeno.substr(pi, 1)      << "\t" 
                     << pi + 1                      << endl;
                if(leftp != 0)
                {
                    (*allofp["crossover_pb_reads.fa"])  << "# " 
                         << (*imitr).first              << "\t"
                         << (*imitr).second.matallele   << "\t"
                         << (*imitr).second.patallele   << "\t"
                         << (*imitr).second.checkallele << "\t" 
                         << thisgeno.substr(pi, 1)      << "\t" 
                         << pi + 1                      << endl;   
                }                 
                pi ++;
                imitr ++;
            }
        }
        // collect informative reads with details
        std::stringstream readss;
        readss.str("");
        //if(alignedtogrp.size() == 0) cout << "   Warning: this grp name is null: check " << pbname << endl;
        readss << alignedtogrp << ":" << spansecond-firstRefsMatch+1;      
        if(pbreadnames_detail.find(pbname) == pbreadnames_detail.end())
        {
            vector<string> alignvec;
            alignvec.push_back(readss.str());              
            pbreadnames_detail.insert(std::pair<string, vector<string> >(pbname, alignvec));
        }else
        {
            pbreadnames_detail[pbname].push_back(readss.str());
        }        
    }
    samifp.close();
    //
    cout << endl;
    if(numNoCG > 0)
    {
        cout << "   Warning: there are a1="
             << numNoCG 
             << " alignments, totaling v1=" 
             << lenNoCG/1000000000.0 
             << " Gb "
             << " without explicit CIAGR info -- collected in unmapped file."   
             << endl;     
    }
    if(numNoMA > 0)
    {
        cout << "   Warning: there are a2="
             << numNoMA 
             << " alignments being secondary/supplementary alignment, skipped. "
             << endl;
    }
    if(numNoSEQ > 0)
    {
        cout << "   Warning: there are a3="
             << numNoSEQ 
             << " alignments without seq field - secondary/supplementary alignment"
             << " or not passing filters, skipped." 
             << endl;
    }
    cout << "   Info: in total a4=" 
         << pbreadnames.size() 
         << " reads from all=" 
         << numraw                  
         << " aligment lines collected (<-header line not counted; raw alignment including secondary etc - skipped in read sep). " 
         << endl;
    // no lg
    cout << "   Info: number of pb alignment seqs WITHOUT linkage info: "   
         << nonlgNo       
         << ", totaling v2=" 
         << lennonlg/1000000000.0      
         << " Gb" 
         << endl;
    cout << "         (u0=" 
         << pbreadnames_nolg.size() 
         << " unique reads in this cluster of no lg or not grouped. )"    
         << endl;
    // lg
    cout << "   Info: number of pb alignment seqs WITH    linkage info a5="   
         << grpedNo 
         << ", totaling v3=" 
         << lengrped/1000000000.0 
         << " Gb"
         << endl;
    cout << "         (u1="<< pbreadnames_lg.size() + pbreadnames_lg_nomkr.size()      
         << " unique reads in this cluster of lg_mkr+lg_nomkr: it is a sum of reads covering or not covering phased markers)"
         << "; among v3 (with a5), "   
         << endl;
    cout << "         a6="        
         << lgnomkr 
         << " covered no phased markers thus cannot determine P/M cluster - "
         << "alignment seq has been put in both P and M cluster, "
         << endl;
    cout << "         taking a portion of v4=" 
         << lennomkr/1000000000.0 
         << " Gb"
         << " (u2="
         << pbreadnames_lg_nomkr.size() 
         << " unique reads in this cluster of lg_nomkr)"                                     
         << endl                  
         << endl;
    cout << "   Note: u1 included all u2. "                                                         << endl;
    cout << "   Note: v1+v2+v3    = total raw pacbio data of (full: 19.929) Gb. "                   << endl;
    cout << "   Note: a1+a4       = total raw pacbio read number (full: 2,587,118)"                 << endl;  
    cout << "   Note: a1+a2+a3+a4 = all raw alignment num; "                                        << endl
         << "         a4 only gives unique readname numbers; 1 readname may have >=2 alignments.)"  << endl;      
    cout << "   Note: u0+u1       = a4. "                                                           << endl;
    // check common between pbreadnames_lg_nomkr and pbreadnames_lg
    map<string, int>::iterator checkitr;
    map<string, int>::iterator checkitr_end;
    checkitr     = pbreadnames_lg_nomkr.begin();
    checkitr_end = pbreadnames_lg_nomkr.end();
    unsigned long commNum = 0;
    unsigned long commLen = 0;    
    while(checkitr != checkitr_end)
    {
        if(pbreadnames_lg.find((*checkitr).first) != pbreadnames_lg.end())
        {
            commNum += 1;
            commLen += (*checkitr).second;
        }
        checkitr ++;
    }
    cout << "   Info: "
         << commNum
         << " reads common between pbreadnames_lg and pbreadnames_lg_nomkr, totaling in "
         << commLen/1000000000.0 
         << " Gb."
         << endl; 
    //
    // read statistics 
    // map<string, map<string, int> > groupedpbreadnames;    
    cout << endl 
         << "   Info: distribution of pacbio reads in linkage groups: "
         << endl << endl;
    map<string, map<string, int> >::iterator gpbitr;
    map<string, map<string, int> >::iterator gpbitr_end;
    gpbitr     = groupedpbreadnames.begin();
    gpbitr_end = groupedpbreadnames.end();
    while(gpbitr != gpbitr_end)
    {
        cout << "\t" << (*gpbitr).first 
             << "\t" << (*gpbitr).second.size() 
             << "\t";
        map<string, int>::iterator nameitr;
        map<string, int>::iterator nameitr_end;
        nameitr     = (*gpbitr).second.begin();
        nameitr_end = (*gpbitr).second.end(); 
        unsigned long redundcount = 0; // one read may have >1 alignments
        while(nameitr != nameitr_end)
        {
            redundcount += (*nameitr).second;
            nameitr ++;
        }
        cout << redundcount << endl;
        gpbitr ++;
    }    
    // close output files
    map<string, ofstream*>::iterator ofileitr;
    map<string, ofstream*>::iterator ofileitr_end;
    ofileitr     = allofp.begin();
    ofileitr_end = allofp.end();
    while(ofileitr != ofileitr_end)
    {
        ofstream *ofp = (*ofileitr).second;
        (*ofp).close();
        //
        ofileitr ++;
    }   
    // for check purpose
    string ofilenamecheck = tmpfolder + "/zcheck_read_group_distribution.txt";
    ofstream checkofp;
    checkofp.open(ofilenamecheck.c_str(), ios::out);
    if(!checkofp.good())
    {
        cout << "   Error: cannot open file to write read-group info file " << ofilenamecheck << endl;
        return 1;
    }
    checkofp << "#pbreadname\tgroup_assigned:score[;...]\tnum_assignment\tbest_assigned_group" << endl;
    map<string, vector<string> >::iterator rgitr;
    map<string, vector<string> >::iterator rgitr_end;
    rgitr     = pbreadnames_detail.begin();
    rgitr_end = pbreadnames_detail.end();
    while(rgitr != rgitr_end)
    {
        unsigned long best_align_score = 0;
        string        best_align_group = "";        
        //cout << "   checking " << (*rgitr).first << endl;        
        checkofp << (*rgitr).first << "\t";
        vector<string> alignvec = (*rgitr).second;
        vector<string>::iterator alitr;
        vector<string>::iterator alitr_end;
        alitr     = alignvec.begin();
        alitr_end = alignvec.end();
        while(alitr != alitr_end)
        {
            checkofp << *alitr;
            //cout << "   checking " << *alitr << endl;
            vector<string> alinfo = split_string(*alitr, ':');
            unsigned long  alisco = strtoul(alinfo[1].c_str(), NULL, 0);
            if(alisco > best_align_score)
            {
                best_align_score = alisco;
                best_align_group = alinfo[0];
            }
            alitr ++;
            if(alitr != alitr_end)
            {
                checkofp << ";";
            }
        }
        checkofp << "\t" << alignvec.size() << "\t" << best_align_group << endl;
        rgitr ++;
    }
    checkofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool overlap_read_del_like(map<unsigned long, DELMARKER>  this_del_like, 
                           unsigned long                  read_sta, 
                           unsigned long                  read_end,
                           map<unsigned long, DELMARKER>* this_del_like_overlapped,
                           string*                        del_phase)
{
    /*
        Input: 
           this_del_like   lists del-like markers along a contig
           read_sta        current pacbio read alignment starting position along the reference 
           read_end        current pacbio read alignment ending   position along the reference            
        Output:            del-like markers overlapping with the current pacbio read.    
    */
    map<unsigned long, DELMARKER>::iterator delitr;
    map<unsigned long, DELMARKER>::iterator delitr_end;
    delitr     = this_del_like.begin();
    delitr_end = this_del_like.end();
    *del_phase = "";
    // score is the overlapping length between del-like and read
    // if both genotypes are collected, we choose the one with larger overlapping region.
    map<string, unsigned long> genotype_max_overlap; 
    while(delitr != delitr_end)
    {
        unsigned long del_sta = (*delitr).second.sta;
        unsigned long del_end = (*delitr).second.end;        
        unsigned long olspan = 0; // overlapping span
        if(del_sta<=read_sta && read_sta<=del_end)
        {
            /*
                      sta---------del---------------end         del-like region
                            sta-----read-----end                pacbio read case 1: right_min = read_end
                            sta-----read------------------end   pacbio read case 2: right_min = del_end
            */
            unsigned long right_min = read_end<del_end?read_end:del_end;
            olspan = right_min - read_sta + 1;
        }else
        if(del_sta<=read_end && read_end<=del_end)
        {
            /*
                      sta---------del---------------end         del-like region
                sta--------------read-----end                   pacbio read case 2: right_min = del_end
            */              
            unsigned long left_max = read_sta>del_sta?read_sta:del_sta;
            olspan = read_end - left_max + 1;
        }else
        if(del_sta>=read_sta && read_end>=del_end)
        {
            /*
                      sta---------del---------------end         del-like region
                sta--------------read---------------------end   pacbio read case 2: right_min = del_end
            */    
            olspan = del_end - del_sta + 1;         
        }
        // do we need a cutoff for considering overlapping as sufficient? current as least 500 bp.
        if(olspan >= 500 || (olspan>=250 && read_end - read_sta + 1 <= 1000))
        {
            (*this_del_like_overlapped).insert(std::pair<unsigned long, DELMARKER>((*delitr).first, (*delitr).second));
            cout << "          : del-like " 
                 << (*delitr).second.ctg    << ":" 
                 << del_sta                 << "-"
                 << del_end                 << " as "
                 << (*delitr).second.geno   << " "
                 << (*delitr).second.haphom << " overlap pb read spanning "
                 << read_sta                << "-" 
                 << read_end                << " by " 
                 << olspan                  << " bp. "
                 << endl;                 
            if((*delitr).second.haphom.find("hap") != std::string::npos)
            {
                *del_phase += (*delitr).second.geno;
                if(genotype_max_overlap.find((*delitr).second.geno) == genotype_max_overlap.end())
                {
                    genotype_max_overlap.insert(std::pair<string, unsigned long>((*delitr).second.geno, olspan));
                }else
                {
                    genotype_max_overlap[(*delitr).second.geno] += olspan;
                }
            }
        }else
        {
            if(0)
            cout << "          : del-like " 
                 << (*delitr).second.ctg    << ":" 
                 << del_sta                 << "-" 
                 << del_end                 << " as "
                 << (*delitr).second.geno   << " "                 
                 << (*delitr).second.haphom << " does not overlap pb read spaning "                               
                 << read_sta                << "-" 
                 << read_end                << "." 
                 << endl; 
        }
        delitr ++;
    }
    //
    if((*del_phase).find("M") != std::string::npos &&
       (*del_phase).find("P") != std::string::npos)
    {
        cout << "          : ambiguious more phasing states found: "            << *del_phase     << endl;
        cout << "            M: " << genotype_max_overlap["M"]                  << " bp versus P: " 
             <<                      genotype_max_overlap["P"]                  << " bp "
             << " (length might have been accumulated from hap-like regions): " << endl;
        if(genotype_max_overlap["M"] > genotype_max_overlap["P"])
        {
            *del_phase = "M";
            cout << "          : however, read would go to genotype " << *del_phase << " with larger overlap. " << endl;
        }else
        if(genotype_max_overlap["M"] < genotype_max_overlap["P"])
        {
            *del_phase = "P";
            cout << "          : however, read would go to genotype " << *del_phase << " with larger overlap. " << endl;
        }else
        {
            *del_phase = "U";
            cout << "          : however, read would go to either genotype P or M due to equal overlap. " << endl;
        }
    }else
    if((*del_phase).find("M") != std::string::npos)
    {
        cout << "          : read would go to genotype " << *del_phase << "."  << endl;
    }else
    if((*del_phase).find("P") != std::string::npos)
    {
        cout << "          : read would go to genotype " << *del_phase << ". " << endl;    
    }else 
    {
        *del_phase = "U";      
        cout << "          : read would go to genotype " << *del_phase << " (either genotype P or M)." << endl;  
    }
    //
    return true;
}
//
bool genotype_pbread(string*                   updated_read, 
                     map<unsigned long, PAMA>* this_marker, 
                     unsigned long             start, 
                     unsigned long             end, 
                     string*                   geno,
                     map<unsigned long, PAMA>* intersectMkr)
{
    (*geno) = "";
    map<unsigned long, PAMA>::iterator positr;
    map<unsigned long, PAMA>::iterator positr_end;
    positr     = (*this_marker).begin();
    positr_end = (*this_marker).end();
    while(positr != positr_end)
    {
        unsigned long thispos = (*positr).first;
        if(thispos>=start && thispos<=end)
        {
            unsigned long this_dist           = thispos - start;
            unsigned long read_correspondence = 0 + this_dist;
            string        this_read_base      = (*updated_read).substr(read_correspondence, 1);            
            string        this_mat_base       = (*positr).second.matallele;
            string        this_pat_base       = (*positr).second.patallele;
            (*positr).second.checkallele      = this_read_base; // put some more info on marker struct
            //
            if(this_read_base.compare(this_mat_base) == 0)
            {
                (*geno) += "M";
            }else
            if(this_read_base.compare(this_pat_base) == 0)
            {
                (*geno) += "P";
            }else
            {
                (*geno) += "u";
            }  
            (*intersectMkr).insert(std::pair<unsigned long, PAMA>((*positr).first, (*positr).second));                           
        }else
        if(thispos>end) 
        {
            break;
        }
        else ;
        //
        positr ++;
    }
    //    
    return true;
}
//
bool align_read_to_ref(vector<char> operation, vector<int> count, string* orignal_read, string* updated_read)
{
    /*
        align the reads to reference according to cigar info:
        //          firstRefsMatch|
        //                        x
        //                        |4M  3D  6M  2I    11M    2D 4M  2d 3M 3S                           CIGAR
        // GCTATTTTGCGGAACAACGAATTCCTGGATCCACCA--GAAAATAAAGTTTGCTGCAGGACTTTTGCCAGCCATAAGCTTCGGTCAGGCT REF
        //                        CCTG---CCACCAGAGAAAATGAAGT--GTTGC--GACT                             READ
        //                        ||||DDD||||||II|||||||||||DD|R|||DD||||                             MM/INDELs
        //                        x         |                        |  x
        //          firstReadMatch|                                     |(secondReadMatch)

        The read would become: 
        //                        CCTG---CCACCA  GAAAATGAAGT--GTTGC--GACT                             READ        
                                               GA
        //
	https://samtools.github.io/hts-specs/SAMv1.pdf: 
        //
	Op	BAM	Description						Consumes_query	Consumes_reference	
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes	
	I	1	insertion to the reference				yes		no	
	D	2	deletion from the reference				no		yes	
	N	3	skipped region from the reference			no		yes	
	S	4	soft clipping (clipped sequences present in	SEQ)	yes		no	
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	no		no	
	P	6	padding (silent deletion from padded reference)		no		no	
	=	7	sequence match						yes		yes	
	X	8	sequence mismatch					yes		yes  	                                                       
    */
    (*updated_read) = "";         // delete 'SI', add 'MDX=' on real read; no action with 'H'; 'NP'
    unsigned long next_start = 0; // on real read
    for(int oi=0; oi<operation.size(); oi++)
    {       
        if(operation[oi]=='H') 
        {
            // need not change read sequence
            (*updated_read) += "";
            next_start      += 0;
        }else
        if(operation[oi]=='S') 
        {
            // clip the read sequence 
            (*updated_read) += "";
            next_start      += count[oi];
        }else
        if(operation[oi]=='M' || operation[oi]=='X' || operation[oi]=='=') 
        {
            (*updated_read) += (*orignal_read).substr(next_start, count[oi]);
            next_start      += count[oi];
        }else
        if(operation[oi]=='D') 
        {
            // extend the read with "-"
            std::string s(count[oi], '-');            
            (*updated_read) += s;
            next_start      += 0;
        }else
        if(operation[oi]=='I') 
        {
            // delete the read subsequence - so update nothing
            (*updated_read) += "";
            next_start      += count[oi];
        }else
        {
            ; // N and P
        }            
    }    
    return true;
}
// read markers with the control by linkage group phasing
bool read_asPollinator_phased_marker(string fmarker,
                                     map<string, bool> ctgJoinMapPhasing_all, 
                                     map<string, map<unsigned long, PAMA> >* correctedmkr)
{
    /*
      format: corr contig_id          position maternal paternal pre-maternal pre-paternal
              corr tig00000013_pilon  21498    C        G        C            G
    */
    ifstream ifp;
    ifp.open(fmarker.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open marker file " << fmarker << endl;
        return false;
    }
    map<string, bool>::iterator cpitr;
    map<string, map<unsigned long, PAMA> >::iterator cmitr;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]  =='#') continue;        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5) 
        {
            cout << "   Warning: unexpected line: " << line;
            continue;
        }
        string        thisctg   = lineinfo[1];
        unsigned long thispos   = strtoul(lineinfo[2].c_str(), NULL, 0);
        string        matallele = lineinfo[3]; // M
        string        patallele = lineinfo[4]; // P
        //
        cpitr = ctgJoinMapPhasing_all.find(thisctg);
        if(cpitr != ctgJoinMapPhasing_all.end())
        {
            if( (*cpitr).second == true)
            {
                string tmpallele = matallele;
                matallele        = patallele;
                patallele        = tmpallele;
            }
        }
        PAMA tpama;
        tpama.matallele = matallele;
        tpama.patallele = patallele;
        //
        cmitr = (*correctedmkr).find(thisctg);
        if(cmitr != (*correctedmkr).end())
        {
            (*cmitr).second.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
        }else
        {
            map<unsigned long, PAMA> posmkr;
            posmkr.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
            (*correctedmkr).insert(std::pair<string, map<unsigned long, PAMA> >(thisctg, posmkr));
        }
    }
    ifp.close();
    return true;
}
//
bool get_del_like_phase(string                                fmarkerdel, 
                 map<string, map<unsigned long, DELMARKER> >* del_like_marker,                      
                 map<string, string>*                         ctgGroup)
{
    /*
        del_like_marker		<contig_id, <start, {ctg, sta, end, haphom, geno, gphase, gmap}> >
        ctgGroup		<contig_id, this_linksage_file>
    */
    //
    cout << "   Info: before reading del-like markers, " 
         << (*ctgGroup).size() << " contigs with known phased grouping. " 
         << endl;
    //
    fstream ifp;
    ifp.open(fmarkerdel.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open del-like marker file " << fmarkerdel << endl;
        return false;
    }
    unsigned long delnum = 0;
    string specialcontig = "";
    map<string, int> unphasedcontig;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        /*
            tig00411999_pilon   
            63737   
            116002  
            hap_low 
            genotype:P   
            phasing_score:3096    
            grouping:8.txt   
            genetic_map:map_group8.txt
        */
        vector<string> lineinfo = split_string(line, '\t');
        DELMARKER tmpdel;
        tmpdel.ctg   = lineinfo[0];        
        tmpdel.sta   = strtoul(lineinfo[1].c_str(), NULL, 0);
        tmpdel.end   = strtoul(lineinfo[2].c_str(), NULL, 0);        
        tmpdel.haphom= lineinfo[3];
        tmpdel.geno  = lineinfo[4];  
        tmpdel.gphase= lineinfo[6];
        tmpdel.gmap  = lineinfo[7];
        // caution: special cases
        // update special cases with contig id changed in marker list
        string thisctg = tmpdel.ctg;
        if(thisctg.find("tig00003826_pilon")!=std::string::npos)
        {
            unsigned long separatepos = 220000; // 210579
            if(tmpdel.end < separatepos)
            {
                thisctg = "tig00003826_pilonsA";                 
            }else
            if(tmpdel.sta > separatepos)
            {
                thisctg = "tig00003826_pilonsB";
            }else
            if(separatepos-tmpdel.sta >= tmpdel.end-separatepos)
            {
                thisctg = "tig00003826_pilonsA";
            }else
            {
                thisctg = "tig00003826_pilonsB";
            }
            if(specialcontig.find(thisctg) == std::string::npos)
            {            
                cout << "          Warning: special contig " << thisctg << endl;
                specialcontig += "\t";
                specialcontig += thisctg;
            }
        }else
        if(thisctg.find("tig00003688_pilon")!=std::string::npos)
        {
            unsigned long separatepos=2050000; 
            if(tmpdel.end < separatepos)
            {
                thisctg = "tig00003688_pilonsA";                 
            }else
            if(tmpdel.sta > separatepos)
            {
                thisctg = "tig00003688_pilonsB";
            }else
            if(separatepos-tmpdel.sta >= tmpdel.end-separatepos)
            {
                thisctg = "tig00003688_pilonsA";
            }else
            {
                thisctg = "tig00003688_pilonsB";
            }
            if(specialcontig.find(thisctg) == std::string::npos)
            {            
                cout << "          Warning: special contig " << thisctg << endl;
                specialcontig += "\t";
                specialcontig += thisctg;
            }
        }else ;   
        // update  
        tmpdel.ctg = thisctg;    
        // collect del-like marker
        if((*del_like_marker).find(tmpdel.ctg) == (*del_like_marker).end())   
        {
            map<unsigned long, DELMARKER> this_del_like;
            this_del_like.insert(std::pair<unsigned long, DELMARKER>(tmpdel.sta, tmpdel));
            (*del_like_marker).insert(std::pair<string, map<unsigned long, DELMARKER> >(tmpdel.ctg, this_del_like));
            delnum ++;
        }else
        {
            (*del_like_marker)[tmpdel.ctg].insert(std::pair<unsigned long, DELMARKER>(tmpdel.sta, tmpdel));
            delnum ++;            
        }
        // update ctg phasing variable
        if( (*ctgGroup).find(tmpdel.ctg) == (*ctgGroup).end() && 
            tmpdel.gphase.find("unphased.txt") == std::string::npos )
        {
            (*ctgGroup).insert(std::pair<string, string>(tmpdel.ctg, tmpdel.gphase));
        }else
        if(tmpdel.gphase.find("unphased.txt") != std::string::npos)
        {
            if(unphasedcontig.find(tmpdel.ctg) == unphasedcontig.end())
            {
                unphasedcontig.insert(std::pair<string, int>(tmpdel.ctg, 1));
            }
        }
        else ;
    }
    vector<string> specialcontiginfo = split_string(specialcontig, '\t');
    if(unphasedcontig.size() > 0)
    {
        cout << "   Info: still, there are " << unphasedcontig.size() << " contigs without grouping/phasing. " << endl;
    }
    cout << "   Info: " << (*ctgGroup).size() 
         << " contigs updated with phasing/grouping according to del-like markers" 
         << " (with " << specialcontiginfo.size()/2 << " contigs become separated) - " 
         << endl;    
    cout << "         originally this is related to " << (*ctgGroup).size() - specialcontiginfo.size()/2 
         << " contigs."
         << endl;  
    cout << "   Info: in total, " << unphasedcontig.size() + (*ctgGroup).size() - specialcontiginfo.size()/2
         << " have been investigated. " 
         << endl;
    /*
    map<string, string>::iterator cpitr;
    map<string, string>::iterator cpitr_end;
    cpitr     = (*ctgGroup).begin();
    cpitr_end = (*ctgGroup).end();
    while(cpitr != cpitr_end)
    {
        cout << "   Check contig group : " << (*cpitr).first << "\t" << (*cpitr).second << endl;
        cpitr ++;
    }
    */
    cout << "   Info: " <<  delnum            << " del-like markers collected. "                  << endl << endl;                        
    //
    ifp.close();
    return true;
}                 
//
bool get_contig_phase(string               jmfile, 
                      map<string, bool>*   ctgJoinMapPhasing,
                      map<string, bool>*   ctgJoinMapPhasing_all,
                      map<string, string>* ctgGroup)
{
    // get ids with the respective phasing status of effective contigs in pollens
    ifstream jmifp;
    jmifp.open(jmfile.c_str(), ios::in);
    if(!jmifp.good())
    {
        cout << "   Error: cannot open file " << jmfile << endl;
        return 1;
    }
    //
    vector<string> this_file_info = split_string(jmfile, '/');
    string jmflagstring = this_file_info[this_file_info.size()-1];
    /*
        414	tig00004185_ri	{0}	...
        413	tig00004185_le	{0}	...
        165	tig00003747_le	{1}	...
    */
    while(jmifp.good())
    {
        string line("");
        getline(jmifp, line);
        if(line.size() == 0 ) continue;
        if(line[0]=='#')      continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<3) 
        {
            cout << "   Warning: unexpected line: " << endl;
            continue;
        }
        //
        if(line.find("{0}")!=std::string::npos || line.find("{1}")!=std::string::npos)
        {
            size_t postig = line.find("tig");
            size_t posund = line.find("_", postig);
            
            string tigid  = line.substr(postig, posund-1-postig+1);
            
            string thiskey = "";
            bool   flip    = true;                        
            if(line.find("le")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon";
                flip    = true;                    
            }else
            if(line.find("le")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon";
                flip    = false;                    
            }else
            if(line.find("ri")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon";
                flip    = true;                    
            }else
            if(line.find("ri")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon";
                flip    = false;                    
            }else ;            
            // caution: special cases
            if(thiskey.find("tigsA003826")!=std::string::npos)
            {
                thiskey = "tig00003826_pilonsA";
                cout << "         special contig: " << thiskey << endl;
            }else
            if(thiskey.find("tig00003826")!=std::string::npos)
            {
                thiskey = "tig00003826_pilonsB";
                cout << "         special contig: " << thiskey << endl;                
            }else            
            if(thiskey.find("tigsA003688")!=std::string::npos)
            {
                thiskey = "tig00003688_pilonsA";
                cout << "         special contig: " << thiskey << endl;                
            }else
            if(thiskey.find("tig00003688")!=std::string::npos)
            {
                thiskey = "tig00003688_pilonsB";
                cout << "         special contig: " << thiskey << endl;                
            }else ;
            // single group-wise
            map<string, bool>::iterator tmpcitr = (*ctgJoinMapPhasing).find(thiskey);
            if(tmpcitr != (*ctgJoinMapPhasing).end())
            {
                if((*tmpcitr).second != flip)
                {
                    cout << "   Warning1: left and right phasings differ on the same contig " << thiskey << endl;
                    cout << "            need to update my code!! "                           << endl;
                    // this should only happens on mis-assembled contigs?
                }
            }else
            {
                (*ctgJoinMapPhasing).insert(std::pair<string, bool>(thiskey, flip));
            }
            // all groups-wise
            tmpcitr = (*ctgJoinMapPhasing_all).find(thiskey);
            if(tmpcitr != (*ctgJoinMapPhasing_all).end())
            {
                if((*tmpcitr).second != flip)
                {
                    cout << "   Warning2: phasings differ on the same contig for "  << thiskey << endl;
                    cout << "            need to update my code!! "                            << endl;
                    // this should only happens on mis-assembled contigs?
                }
            }else
            {
                (*ctgJoinMapPhasing_all).insert(std::pair<string, bool>(thiskey, flip));
            }     
            
            map<string, string>::iterator c2jmitr;
            c2jmitr = (*ctgGroup).find(thiskey);
            if(c2jmitr != (*ctgGroup).end())
            {
                if((*c2jmitr).second.compare(jmflagstring) != 0)
                {
                    cout << "   Warning3: same contig found in different linkage groups? It can happen! " << endl;
                    //
                    (*c2jmitr).second += "#";
                    (*c2jmitr).second += jmflagstring; // later need to send related reads to both group.
                }
            }
            else
            {
                (*ctgGroup).insert(std::pair<string, string>(thiskey, jmflagstring));            
            }       
            //
        }
    }
    jmifp.close();
    
    return true;
}
//
bool get_contig_phase_files(string fphase, vector<string>* groupphase)
{
    ifstream ifp;
    ifp.open(fphase.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << fphase << endl;
        return false; 
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]  =='#') continue;
        
        ifstream tmpifp;
        tmpifp.open(line.c_str(), ios::in);
        if(!tmpifp.good())
        {
            cout << "   Warning: group phase file link not valid, skipped: " << line << endl;
            
        }else
        {
            (*groupphase).push_back(line);
            cout << "   Info: group phase file collected: "                  << line << endl;
            tmpifp.close();
        }
    }
    ifp.close();
    if((*groupphase).size()==0) return false;
    else return true;
}
//
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count)
{
    /*
	https://samtools.github.io/hts-specs/SAMv1.pdf: 
        //
	Op	BAM	Description						Consumes_query	Consumes_reference	
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes	
	I	1	insertion to the reference				yes		no	
	D	2	deletion from the reference				no		yes	
	N	3	skipped region from the reference			no		yes	
	S	4	soft clipping (clipped sequences present in	SEQ)	yes		no	
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	no		no	
	P	6	padding (silent deletion from padded reference)		no		no	
	=	7	sequence match						yes		yes	
	X	8	sequence mismatch					yes		yes  	
	NOTE:
        // CIGAR alphabet: [0-9]+[MIDNSHPX=] and *: 23S23M1D8M1I16M1D29M38H 		
        // Sum of lengths of the MIS=X operations shall equal the length of SEQ in bam 
        // H can only be present as the first and/or last operation.
        // S may only have H operations between them and the ends of the CIGAR string
        // or mRNA-to-genome alignment, an N operation represents an intron.  	  
    */
    char *cstr = (char*)cigar.c_str(); // 
    string numstr("");
    int covlen = 0;
    for (int i=0; i<cigar.size(); i++)
    {
        if(cstr[i]>='0' && cstr[i]<='9')
        {
            numstr += cstr[i];
        }
        else
        if(cstr[i] == 'M' || cstr[i] == 'I' || cstr[i] == 'D' || 
           cstr[i] == 'N' || 
           cstr[i] == 'S' || cstr[i] == 'H' || 
           cstr[i] == 'P' || 
           cstr[i] == 'X' || cstr[i] == '=')
        {
            (*operation).push_back(cstr[i]);    
            (*count).push_back(atoi(numstr.c_str()));
            //                    
            if(cstr[i] == 'M' || cstr[i] == 'D' || cstr[i] == 'N' || cstr[i] == 'X' || cstr[i] == '=')
            {
                covlen += atoi(numstr.c_str()); // positions of ref covered by read
            }
            numstr.clear();
        }
        if(cstr[i] == 'X' || cstr[i] == '=' || cstr[i] == 'P' || cstr[i] == 'N') 
        {            
            cout << "   Warning: you have special operation \'" << cstr[i] << "\' in Cigar: " << cigar << endl;
        }
    }
    // check
    if(true)
    {
        //cout << endl << "   check: cigar=" << cigar << endl;
        for(int ci = 0; ci < (*operation).size(); ci ++)
        {
            if(ci!=0 && ci!=(*operation).size()-1)
            {
                char tmpc = (*operation)[ci];
                if(tmpc=='S' || tmpc=='H')
                {
                    cout << "   Warning: \'" << tmpc << "\' operation happened in middle of alignment with cigar " 
                         << cigar            << endl;
                }            
            }
            // cout << "   check: operation " << (*operation)[ci] << "\t" << (*count)[ci] << endl;
        }
        // cout << "   check: reference length covered: " << covlen << " bp. " << endl; 
    }
    //
    return covlen;
}
/* 0.QNAME 1.FLAG 2.RNAME 3.POS 4.MAPQ 5.CIGAR 6.RNEXT 7.PNEXT 8.TLEN 9.SEQ 10.QUAL

 0.QNAME:   M05453:196:000000000-CGNKH:1:1101:16586:1168    
 1.FLAG:    163 
 2.REFNAME: refname    
 3.POS:     1   
 4.MAPQ:    42  
 5.CIGAR:   129M    
 6.RNEXT:   =   
 7.PNEXT:   208 
 8.REFLEN:  431 
 9.SEQ:     GGTCAAGGCAAGACGATATAACTGAACTCCGTTGTAGCATTAGAGCTGAAATGTTCTGTGGTTGAATTAATTTGTTTCTGGCAAATAATTAAAGTTGTTGCTGTTGGATTTACGTTGTAGGTATTTGGG   
10.QUAL:    GF9@EFFFCFGGGGG+,CFCCCFFFFFGCGED@GGCCF<ACE96CECG@<,CFF<FF@6EFC8FEGDFGFCFFCC@,CE@EE@F8EF9FFC<@,CEEFF8FFF:B<8AFFGGG,,BEGFGGFEAFFC?7   
 ...
*/
//
//
bool get_options(int     argc,
                 char*   argv[],
                 string* fsam,
                 string* fmarker,
                 string* fmarkerdel,
                 string* fphase,
                 string* outprefix,
                 bool*   verbose,
                 double* minCOscore)
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
        if(optstr.compare("--sam") == 0)
        {
            ic ++;
            *fsam = (string)argv[ic];
            ifstream fp;
            fp.open((*fsam).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: SAM file provided: "         
                     << *fsam << endl;
            }
            else
            {
                cout << "   Error: cannot open SAM file "      
                     << *fsam << endl;
                return false;
            }
        }
        else
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
            }
            else
            {
                cout << "   Error: cannot open marker file "      
                     << *fmarker << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--marker2") == 0)
        {
            ic ++;
            *fmarkerdel = (string)argv[ic];
            ifstream fp;
            fp.open((*fmarkerdel).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: del-like marker file provided: "         
                     << *fmarkerdel << endl;
            }
            else
            {
                cout << "   Error: cannot open del-like marker file "      
                     << *fmarkerdel << endl;
                return false;
            }
        }        
        else
        if(optstr.compare("--phase") == 0)
        {
            ic ++;
            *fphase = (string)argv[ic];
            ifstream fp;
            fp.open((*fphase).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: contig phasing file provided: "         
                     << *fphase << endl;
            }
            else
            {
                cout << "   Error: cannot open contig phasing file "      
                     << *fphase << endl;
                return false;
            }
        }else
        if(optstr.compare("--ims") == 0)
        {
            ic ++;
            *minCOscore = (double)atof(argv[ic]);
            if(*minCOscore>=0.5)
            {
                cout << "   Info: minimum score of candidate co events provided: "   << *minCOscore << endl;
            }
            else
            {
                cout << "   Error: minimum score of candidate co events too small: " << *minCOscore << endl;
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
//
string get_break_pos(string           poPat, 
                     string*          pmflag,
                     unsigned long*   leftp, 
                     double*          score)
{
    // used after 20191110
    // poPat: is the PMu pattern sequence
    // leftp : would be the i-th   marker position along the PMu sequence
    // score : would be score for the reported co
    // minimum number of markers hard-coded
    int min_marker = 8;
    //
    string pmstring = "";
    string poPatbk  = poPat;
    poPatbk.erase(std::remove(poPatbk.begin(), poPatbk.end(),'u'), poPatbk.end()); // remove non-effective marekers with value: 'u'
    size_t pcnt = std::count(poPatbk.begin(), poPatbk.end(), 'P');
    size_t mcnt = std::count(poPatbk.begin(), poPatbk.end(), 'M');
    if(pcnt > mcnt)
    {
        *pmflag = "P";
    }else
    if(mcnt > pcnt)
    {
        *pmflag = "M";
    }else
    {
        if(mcnt > 0)
        {
            srand (time(NULL));
            double randra = rand()%1000000000*1.0 / 1000000000.0;   
            if(randra <= 0.5)
            {
                *pmflag = "M";
            }else
            {
                *pmflag = "P";
            }
        }else
        {
            *pmflag = "U"; // go to both (indeed, either M or P - randomly assignment, so expected 50% to each)
        }
    }
    // too few markers
    if(poPatbk.size()<min_marker) // do have consider CO with limited markers 2019-11-10
    {
        *leftp   = 0; // no co
        *score   = 0; // no co         
        pmstring = (*pmflag) + (*pmflag);
        cout << "   Check: " << " ---- determined as "  << pmstring << " within initial cluster " << *pmflag 
             << " ( " << poPatbk.size() << " effective snp marker - too few) "                    << endl;         
        return pmstring; 
    }
    // pure P
    if(mcnt==0 && pcnt>0)
    {
        *leftp   = 0; // no co
        *score   = 0; // no co         
        pmstring = (*pmflag) + (*pmflag);
        cout << "   Check: " << " ---- determined as "  << pmstring << " within initial cluster " << *pmflag 
             << " ( " << poPatbk.size() << " effective snp marker - pure) "                       << endl;         
        return pmstring;         
    }
    // pure M
    if(pcnt==0 && mcnt>0)
    {
        *leftp   = 0; // no co
        *score   = 0; // no co         
        pmstring = (*pmflag) + (*pmflag);
        cout << "   Check: " << " ---- determined as "  << pmstring << " within initial cluster " << *pmflag 
             << " ( " << poPatbk.size() << " effective snp marker - pure) "                       << endl;         
        return pmstring;         
    }    
    // PMu for below cases:
    *score        = 0.0;
    size_t c1Lef  = 0;   // Left  P cnt for highest score
    size_t c2Lef  = 0;   // Left  M cnt for highest score
    size_t c1Rig  = 0;   // right P cnt for highest score
    size_t c2Rig  = 0;   // right M cnt for highest score
    for(size_t mi = 1; mi <= poPat.size()-1; mi ++) // later control on number of markers, and no need checking all pos.
    {
        string strLef = poPat.substr(0, mi); // val @ pos of 0,...,mi-1
        string strRig = poPat.substr(mi);    // val @ pos of mi,mi+1,...
        // current break point: note: 'u' excluded from consideration
        size_t c1Leftmp  = std::count(strLef.begin(), strLef.end(), 'P'); // Left  P
        size_t c2Leftmp  = std::count(strLef.begin(), strLef.end(), 'M'); // Left  M
        size_t c1Rigtmp  = std::count(strRig.begin(), strRig.end(), 'P'); // Right P
        size_t c2Rigtmp  = std::count(strRig.begin(), strRig.end(), 'M'); // Right M
        //
        double af1Leftmp = (double)c1Leftmp*1.0/(c1Leftmp+c2Leftmp);
        double af2Leftmp = (double)c2Leftmp*1.0/(c1Leftmp+c2Leftmp);        
        double af1Rigtmp = (double)c1Rigtmp*1.0/(c1Rigtmp+c2Rigtmp);
        double af2Rigtmp = (double)c2Rigtmp*1.0/(c1Rigtmp+c2Rigtmp);
        //
        double score1 = af1Leftmp * af2Rigtmp; // left as P right as M
        double score2 = af2Leftmp * af1Rigtmp; // left as M right as P
        double scotmp = score1>score2?score1:score2;   
        if(scotmp >= *score)
        {
            *score  = scotmp;
            *leftp  = mi;  // (*leftp)-th marker is the left break marker, (*leftp+1) is the right break marker
            if(score1>score2)
            { 
                pmstring  = "PM";
            }else
            {
                pmstring  = "MP";
            }  
            //
            c1Lef  = c1Leftmp;
            c2Lef  = c2Leftmp;
            c1Rig  = c1Rigtmp;
            c2Rig  = c2Rigtmp;                     
        }
    }
    cout << "   Check: best CO candidate after " << *leftp << "th marker. " << endl;                        
    // 
    string checkfrom = "";
    // control on marker number: MuMMMMuuuuuuPuPPPu
    string lowsco = "";
    if(*score>=minCOscore)
    {
        // from cluster "M" 
        if((*pmflag).compare("M")    == 0)
        {        
            // but found "P" on left end 
            if(pmstring.compare("PM") == 0)
            {            
                // and at least 5 effective "P" and first position was "P"
                if((poPatbk.substr(0, 1)).compare("P") == 0 && c1Lef >= min_marker/2 )
                {
                    checkfrom = "9";
                }else
                {
                    pmstring = "MM";
                    *leftp   = 0;  // no co
                    //*score  = 0; // no co
                    checkfrom = "10";
                }
            }
            // but found "P" on right end
            if(pmstring.compare("MP") == 0)
            {
                // and at least 5 effective "P" and last position was "P"
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("P") == 0 && c1Rig >= min_marker/2 )
                {
                    checkfrom = "11";
                }else
                {
                    pmstring = "MM";
                    *leftp   = 0;  // no co
                    //*score  = 0; // no co 
                    checkfrom = "12";
                }
            }
        }
        // from cluster "P"
        if((*pmflag).compare("P")    == 0)
        {        
            // but found "M" on left end 
            if(pmstring.compare("MP") == 0)
            {            
                // and at least 5 effective "M" and first position was "M"                
                if((poPatbk.substr(0, 1)).compare("M") == 0 && c2Lef >= min_marker/2)
                {
                    checkfrom = "13";
                }else
                {
                    pmstring = "PP";
                    *leftp   = 0; // no co
                    //*score  = 0; // no co
                    checkfrom = "14";
                }
            }
            // but found "M" on right end and effective last  position was "M"
            if(pmstring.compare("PM") == 0)
            {            
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("M") == 0 && c2Rig >= min_marker/2)
                {
                    ;
                }else
                {
                    pmstring = "PP";
                    *leftp   = 0; // no co
                    //*score  = 0; // no co
                    checkfrom = "15";
                }
            }
        }       
    }else
    {
        lowsco = "low ";
        if((*pmflag).compare("M") == 0)
        {
            pmstring     = "MM";
            checkfrom = "16";
        }else
        if((*pmflag).compare("P") == 0)
        {
            pmstring     = "PP";
            checkfrom = "17";
        }
        *leftp   = 0; // no co
        //*score = 0; // no co
    }
    //
    if(*leftp   == 0)
    {
        pmstring = (*pmflag) + (*pmflag);
    }
    //
    if(pmstring.compare("PM")==0 || pmstring.compare("MP")==0)
    {
        cout << "   Check: " << " ---- determined as "  << pmstring << " within given cluster "    << *pmflag 
             << " (CO after "<< *leftp << "th marker: " << lowsco   << "score "                    << *score   
             << " - case "   << checkfrom << " ) "      << endl;
    }else
    {
        cout << "   Check: "    << " ---- determined as "  << pmstring << " within given cluster " << *pmflag 
             << " (no CO after "<< *leftp << "th marker: " << lowsco   << "score "                 << *score   
             << " - case "      << checkfrom << " ) "      << endl;     
    }
    return pmstring;            
}
