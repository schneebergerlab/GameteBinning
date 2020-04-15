/* this function, given a SNP marker file and genotype data of multiple pollen samples
   phase SNPs to reconstruct haplotypes.
   
   species targeted: heterozygous
   
   TODO: consider small indels as markers
   
   v5:  targeted checking like tig00004000_pilon, why there are enough markers but not correctly phased!
   v6:  output crossovers
   v6.1 output crossovers with detailed marker positions
   
   Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de 
   2019-2020
*/ 
#include                   <stdio.h>
#include                  <string.h>
#include                    <string>
#include                       <map>
#include                  <stdlib.h>
#include                  <iostream>
#include                   <sstream>
#include                   <fstream>
#include                 <algorithm>
#include                    <time.h>
#include                  <assert.h>
#include                    <math.h>
#include                   <iomanip>
#include                <sys/stat.h>
#include                  <dirent.h>
#include            "split_string.h"
//
struct ALLELE
{
    string refb;    // given ref might be from A or B
    string altb;    // given ref might be from B or A
    string pollenb; // what is the allele in pollen?
};
struct POLLEN
{
    string poSeq;   // seq at markers for pollen haplotype
    string poAll;   // over all, can be determined as 'M', 'P', 'MP', 'PM', or 'N'
    string poPat;   // string like MMMMMMppMMMMMMMpMMMMMMMM: indicating parental genotypes
    vector<unsigned long> mkrpos; // position of markers for this contig
};
struct NODE
{
    string contig; // current contig x
    NODE    *next; // best contig to the right of contig x
};
struct CLUSCNT
{
    int Pcnt;
    int Mcnt;
    int ucnt;
};
struct CTGLG
{
    double lod[2];
    string edge1;
    string edge2;
};
struct BIPAT
{
    string leftPat;
    string rightPat;
};
double minCOscore  = 0.64; // 
//
bool collect_marker_info(const char*                        file, 
                         map<string, map<unsigned long, ALLELE> >* mkr);
bool collect_contig_size(string                             sizefile, 
                         map<string, unsigned long>*        contigsize);                         
bool get_options(        int                                argc, 
                         char*                              argv[],
                         string*                            filemark,
                         string*                            filesize,
                         vector<string>*                    filepoll,
                         bool*                              correction,
                         double*                            lodcutoff,
                         double*                            minCOscore,
                         string*                            outprefix);
bool get_base_position(  string                             base,                         
                         int*                               pos);
bool get_pollen_allele(  string                             pollenfile, 
                         map<string, map<unsigned long, ALLELE> >* mkr);                         
bool cluster_pollens(    map<int, POLLEN>                   pollens, 
                         map<int, POLLEN>*                  maCluster,
                         map<int, POLLEN>*                  paCluster,
                         map<int, POLLEN>*                  unCluster,
                         map<int, POLLEN>*                  pollens_updated);
bool flip_pcluster(      map<int, POLLEN>*                  paCluster,
                         string                             matGT,
                         string                             patGT);                         
bool get_similarty(      string                             matGT, 
                         string                             patGT, 
                         string                             pollenGT,
                         string*                            pattern,                   
                         unsigned long*                     matscore, 
                         unsigned long*                     patscore);
// correct raw haplotypes from assembly with pollen haplotypes                         
bool correct_haplotype(  map<int, POLLEN>                   pollens,
                         map<int, POLLEN>*                  pollens_corrected);
bool correct_haplotype2_with_consensus(
                         map<int, POLLEN>                   maCluster,
                         map<int, POLLEN>                   paCluster,
                         string*                            matGT,
                         string*                            patGT);
bool correct_haplotype3_with_pmClusterCnt(
                         map<int, POLLEN>                   maCluster,
                         map<int, POLLEN>                   paCluster,
                         string*                            matGT,
                         string*                            patGT);    
bool get_PM_count(       string                             posPattern, 
                         CLUSCNT*                           posCNT);           
double calculate_lod(    string                             patA, 
                         string                             patB);                                                            
bool order_linked_contigs(map<string, BIPAT>                contigPollenGTSeq, 
                         string                             tmpfolder);                         
bool reverse_contig_GT(  string                             contig_GT, 
                         string*                            contig_GT_rev);       
bool get_match_score(    string                             contig_GT1, 
                         string                             contig_GT2, 
                         int*                               score);
string get_crossover(    string                             poPat, 
                         string                             pmflag);
string get_break_pos(    string                             poPat, 
                         string                             pmflag,
                         int                                pollenid,
                         unsigned long*                     leftp, 
                         double*                            score,
                         string                             contigname,
                         map<unsigned long, ALLELE>         contigmarker);   
bool find_bed(           string                             poPat, 
                         string                             pmflag, 
                         int                                pollenid, 
                         bool*                              out, 
                         vector<string>*                    bedregions);
unsigned long recover_real_marker_pos(int                   markerposition,
                         map<unsigned long, ALLELE>         contigmarker);
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 9)
    {
        cout << endl;
        cout << "   Given raw (snp) marker files defined along contigs/chrs and \n" 
             << "         resequencings of multiple pollen genomes, \n"            
             << "   this tool phases markers (and order haplotye contigs - TODO)."           << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "   (version 1.0 - compiled on " << buildString << ")"               << endl << endl;
        cout << "   Usage: asPollinator --marker snp_list.txt --pollen pollen_consensus_flist.txt --size chrsizes.txt [--corr] [--lod] [--ims] -o outprefix" << endl;
        cout << "          --marker gives snp marker file. "                                 << endl;
        cout << "          --pollen gives lists of A/C/G/T count info files. "               << endl;
        cout << "          --size   gives sizes of contigs. "                                << endl;
        cout << "          --corr   gives whether phasing needed. "                          << endl;
        cout << "          --ims    gives the minimum score for defining a co [0.64] "       << endl;
        cout << "          --lod    gives the lod value to build contig connections - TODO." << endl;
        cout << endl;
        exit(1);
    }
    double startT= clock();
    // set variables from options
    // marker file: format: projn chr pos ref alt ...
    string filemark("");
    // contig size file: for creating bed info
    string filesize("");
    // list of pollen files: each line is one path/to/pollenx/quality_variant.txt
    vector<string> filepoll; 
    //
    bool  correction = false; // default no haplotype correction
    double lodcutoff = 40.0;
    // output label
    string outprefix;
    if(!get_options(argc, 
                    argv,
                    &filemark,
                    &filesize,
                    &filepoll,
                    &correction,
                    &lodcutoff,
                    &minCOscore,
                    &outprefix)
      )
    {
        cout << "   Error: incorrect parameter settings. " << endl;
        return false;
    }    
    // step 0. initialized a variable for recording known marker info:  <contig, <pos, {ref, alt}> >
    map<string, map<unsigned long, ALLELE> >  knownMarker;
    // step 1. read known marker info: only single
    if(!collect_marker_info(filemark.c_str(), &knownMarker))
    {
        cout << "   Error: cannot open marker file " << filemark << "; exited." << endl;
        exit(1);
    }
    // read size info
    map<string, unsigned long> contigsize;
    if(filesize.size()==0)
    {
        cout << "   Erorr: contig size file not collected, check your options..." << endl;
        exit(1); 
    }
    if(!collect_contig_size(filesize, &contigsize))
    {
        cout << "   Error: failed in collecting contig size info." << endl;
        exit(1);
    }
    // step 2: get pollen allele at marker: can be 
    //         refb; altb; none of refb and altb found 'N'; both refb and altb found 'U'
    // seq at markers for ctg of pollen sample: <contigID, <pollenID, {Seq, overallGT} > >
    // seq at markers for ctg of maternal and paternal as well; key as "haplotype"
    // create an intermediate folder for collecting details
    string tmpfolder = outprefix+"_tmp_pollen_genotypes";
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
    string tmpfolder_s1 = outprefix+"_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise";
    dir = opendir(tmpfolder_s1.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder_s1.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder_s1 << endl;
            return false;
        }
    }
    else;
    // subfolder for collecting bed information of P/M of pollens along a contig 
    string tmpfolder_s6 = outprefix+"_tmp_pollen_genotypes/s6_PM_pollen_bed_ctgwise";
    dir = opendir(tmpfolder_s6.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder_s6.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder_s6 << endl;
            return false;
        }
    }
    else;
    // file for collecting phased markers
    std::stringstream phasemarkerinfo;
    phasemarkerinfo.str("");
    phasemarkerinfo << tmpfolder
                    << "/s4_phased_markers.txt"; // contigs::both-haplotypes,ith-pollenSeq-at-markers
    ofstream phasedmkrofp;
    phasedmkrofp.open((phasemarkerinfo.str()).c_str(), ios::out);
    if(!phasedmkrofp.good())
    {
        cout   << "   Error: cannot open file " << phasemarkerinfo.str() << endl;
        return false;
    }
    // file for collecting statistics of markers
    string markerStatFile = tmpfolder + "/" + "s5_ctg_markers_stat.txt";    
    ofstream statmkrofp;
    statmkrofp.open(markerStatFile.c_str(), ios::out);    
    if(!statmkrofp.good())
    {
        cout << "   Error: cannot open file " << markerStatFile << endl;
        return false;
    }  
    statmkrofp << "#Ctg-ID\tMarker-Num\tPM-Cluster\tPollen-Num\tTotal-Covered-Marker-Num\tTotal-Confusing-Marker-Num\tPer-Marker-Coverage" 
               << endl;      
    //  
    map<string, map<int, POLLEN> >  pollenCtgSeq;
    //
    vector<string>::iterator pfitr;
    vector<string>::iterator pfitr_end;
    pfitr     = filepoll.begin();
    pfitr_end = filepoll.end();
    int fid   = 0;
    while(pfitr != pfitr_end)
    {
        string this_pofile = *pfitr;
        vector<string> this_pofile_info = split_string(this_pofile, '/');
        fid ++;
        std::stringstream fidstr;
        fidstr.str("");
        fidstr << "po" << fid;     
        //   
        cout << "   Info: reading genotype info for " << fid << "th pollen file..." << endl;   
        if(!get_pollen_allele(this_pofile, 
                              &knownMarker))
        {
            return false;
        }
        // collect genotype-sequences at markers for each contig for each pollen
        std::stringstream gtinfo;
        gtinfo.str("");
        gtinfo << tmpfolder
               << "/genotype_pollen_"
               << fid
               << ".txt\0"; // ordered fid according to inputs; caution on matching files!
        ofstream tmpofp;
        tmpofp.open((gtinfo.str()).c_str(), ios::out);
        if(!tmpofp.good())
        {
            cout   << "   Error: cannot open file " << gtinfo.str() << endl;
            return false;
        }
        //
        map<string, map<unsigned long, ALLELE> >::iterator ctgitr;
        map<string, map<unsigned long, ALLELE> >::iterator ctgitr_end;
        ctgitr     = knownMarker.begin();
        ctgitr_end = knownMarker.end();
        int effective_num = 0;
        while(ctgitr != ctgitr_end)
        {
            string this_ctg_id  = (*ctgitr).first;
            POLLEN tmppollen;
            tmppollen.poSeq = "";
            tmppollen.poAll = "";
            POLLEN  tmphaplo;
            tmphaplo.poSeq  = ""; // for maternal seq
            tmphaplo.poAll  = ""; // for paternal seq
            //
            map<unsigned long, ALLELE>::iterator positr;
            map<unsigned long, ALLELE>::iterator positr_end;
            positr     = (*ctgitr).second.begin();
            positr_end = (*ctgitr).second.end();
            while(positr != positr_end)
            {
                string pgt       = (*positr).second.pollenb;
                tmppollen.poSeq += pgt;
                tmpofp << (*ctgitr).first          << "\t"
                       << (*positr).first          << "\t"
                       << (*positr).second.refb    << "\t"
                       << (*positr).second.altb    << "\t"
                       << pgt                      << endl;
                if(pgt.compare("N")!=0 && pgt.compare("U")!=0) effective_num ++;
                tmphaplo.poSeq += (*positr).second.refb; // initial maternal haplotype
                tmphaplo.poAll += (*positr).second.altb; // initial paternal haplotype
                // marker positions
                // unsigned long mkrpostmp = strtoul( ((*positr).first).c_str(), NULL, 0);
                unsigned long mkrpostmp = (*positr).first;
                (tmphaplo.mkrpos).push_back(mkrpostmp);
                // reset
                (*positr).second.pollenb = "N";
                positr ++;
            }
            // collect seq of this contig for current pollen
            map<string, map<int, POLLEN> >::iterator poctgitr = pollenCtgSeq.find(this_ctg_id);
            if(poctgitr == pollenCtgSeq.end())
            {
                // pollen
                map<int, POLLEN> tmpPINFO;
                tmpPINFO.insert(std::pair<int, POLLEN>(fid, tmppollen)); // <pid, {pseq, pall} >
                pollenCtgSeq.insert(std::pair<string, map<int, POLLEN> >(this_ctg_id, tmpPINFO));
                // maternal and paternal
                poctgitr = pollenCtgSeq.find(this_ctg_id);
                (*poctgitr).second.insert(std::pair<int, POLLEN>(0, tmphaplo));      
            }else
            {
                // pollen
                (*poctgitr).second.insert(std::pair<int, POLLEN>(fid, tmppollen));
                // maternal and paternal -- no need repeating
            }
            // next contig
            ctgitr ++;
        }
        tmpofp.close();
        cout << "      with " << effective_num << " effective markers with clear genotype. " <<  endl;
        //
        pfitr ++;
    }
    //   
    // step 1: checking each contig marker-sequences info
    // for each contig, define a sequence of 'P/M/U' for pollen 1-N; P:maternal haplotype, M: maternal haplotype
    map<string, map<int, string> > contigsPM; // <contigid, <pollenid, pollenContigPM > >
    //  
    std::stringstream gtinfo;
    std::stringstream bedinfo;    
    map<string, map<int, POLLEN> >::iterator pcitr;
    map<string, map<int, POLLEN> >::iterator pcitr_end;
    pcitr     = pollenCtgSeq.begin();
    pcitr_end = pollenCtgSeq.end();
    while(pcitr != pcitr_end)
    {
        //
        map<string, map<unsigned long, ALLELE> >::iterator ccitr;
        ccitr     = knownMarker.find( (*pcitr).first );
        if(ccitr == knownMarker.end())
        {
            cout << "   Error: contig " << (*pcitr).first << " not found in marker list?" << endl;
            return false;
        }
        map<unsigned long, ALLELE> contigmarker = (*ccitr).second; // for creating bed file
        //
        gtinfo.str("");
        gtinfo << tmpfolder_s1 
               << "/s1_genotype_pollen_seq_contig_"
               << (*pcitr).first 
               << ".txt\0"; // contigs::both-haplotypes,ith-pollenSeq-at-markers
        ofstream mergeofp;
        mergeofp.open((gtinfo.str()).c_str(), ios::out);
        if(!mergeofp.good())
        {
            cout   << "   Error: cannot open file " << gtinfo.str() << endl;
            return false;
        }
        //
        bedinfo.str("");
        bedinfo << tmpfolder_s6
                << "/s6_PM_region_pollens_at_contig_"
                << (*pcitr).first 
                << ".txt\0"; // contigs
        ofstream bedofp;
        bedofp.open((bedinfo.str()).c_str(), ios::out);
        if(!bedofp.good())
        {
            cout   << "   Error: cannot open file " << bedinfo.str() << endl;
            return false;
        }
        bedofp << "#contig-id"  << "\t" 
               << "start-pos"   << "\t" 
               << "end-pos"     << "\t" 
               << "start-index" << "\t"
               << "end-index"   << "\t"
               << "region-pat"  << "\t"
               << "PM-cluster"  << "\t"
               << "pollen-id_x" << endl;
        bedofp << "# note: the first region will be left-extended to coordinate 1 on contig " << endl
               << "#       the last region will be right-extended to coordinate contig-size " << endl;        
        //
        string ctgid = (*pcitr).first;
        mergeofp << "#>" << ctgid << endl;
        map<int, POLLEN> pollens   = (*pcitr).second;
        string matGT = pollens[0].poSeq; // initially assembled maternal backbone --------------------------------------
        string patGT = pollens[0].poAll; // initially assembled paternal backbone --------------------------------------
        vector<unsigned long> mkrpos = pollens[0].mkrpos;
        // output initial pollens
        map<int, POLLEN>::iterator ppitr;
        map<int, POLLEN>::iterator ppitr_end;
        ppitr     = pollens.begin();
        ppitr_end = pollens.end();
        while(ppitr != ppitr_end)
        {
            POLLEN ptmp = (*ppitr).second;
            if((*ppitr).first == 0)
            {
                mergeofp << "#" << ptmp.poSeq << "\t" << (*ppitr).first << "_m provided/assembled " << endl;
                mergeofp << "#" << ptmp.poAll << "\t" << (*ppitr).first << "_p provided/assembled " << endl;
            }else
            {
                //mergeofp << ptmp.poSeq << "\t" << (*ppitr).first << "_x" << endl;
            }
            ppitr ++;
        }
        //        
        // sub step 1. correct maternal and paternal haplotypes with pollen sample info: first correction
        //             bi-marker phasing
        cout << "   Check: correct haplotypes for contig " << ctgid << endl; 
        map<int, POLLEN> pollens_corrected;
        cout << "   Check: 1st correcting haplotypes (bi-marker phasing)..."   << endl;
        if(!correction)
        {
            cout << "   Info: provided parental haplotypes not corrected as user required." << endl;
            pollens_corrected = pollens;
        }else   
        if(!correct_haplotype(pollens, &pollens_corrected))
        {
            cout << "   Error: 1st correction failed. " << endl;
            return false;
        }
        else ;
        string matGT_corrected = pollens_corrected[0].poSeq; // _corrected maternal backbone ---------------------------
        string patGT_corrected = pollens_corrected[0].poAll; // _corrected paternal backbone ---------------------------                    
        // 
        // sub step 2. output clustered pollens: pollens_corrected: first clustering
        map<int, POLLEN> maCluster;
        map<int, POLLEN> paCluster;
        map<int, POLLEN> unCluster;
        map<int, POLLEN> pollens_updated;
        // cluster pollens into paternal and maternal; and set attribute at each marker for each pollen
        cout << "   Check: 1st clustering... " << endl;
        if(!cluster_pollens(pollens_corrected, 
                            &maCluster,
                            &paCluster,
                            &unCluster,
                            &pollens_updated)
          )              
        {
            return false;
        }        
        // sub step 3. second correction with flipped the whole pollen sequence and pattern
        if(correction)  
        {      
            cout << "   Check: flipping paCluster... " << endl;
        }
        if(correction && !flip_pcluster(&paCluster, 
                          matGT_corrected,
                          patGT_corrected))
        {
            return false;
        }
        if(correction)
        {
            cout << "   Check: 2nd correcting haplotyes (find consensus with 1st clustering)..." << endl;
        }
        if(correction && !correct_haplotype2_with_consensus(maCluster,
                                                            paCluster,
                                                            &matGT_corrected,
                                                            &patGT_corrected) )
        {
            return false;
        }
        // sub step 4. cluster again with consensus haplotypes
        cout << "   Check: 2nd clustering... " << endl;
        maCluster.clear();
        paCluster.clear();
        unCluster.clear();
        pollens_updated.clear();
        // 
        pollens_corrected[0].poSeq = matGT_corrected;
        pollens_corrected[0].poAll = patGT_corrected;        
        if(!cluster_pollens(pollens_corrected, 
                            &maCluster,
                            &paCluster,
                            &unCluster,
                            &pollens_updated)
          )              
        {
            return false;
        }        
        // sub step 5. correct haplotypes according to P/M counts in paternal and maternal clusters
        if(correction)
        {
            cout << "   Check: 3rd correcting haplotypes (with P/M counts of 2nd clustering)... " 
                 << endl;
        }
        if(correction && !correct_haplotype3_with_pmClusterCnt(maCluster,
                                                               paCluster,
                                                               &matGT_corrected,
                                                               &patGT_corrected) )
        {
            return false;
        }
        //        
        // sub step 6. 3rd cluster with final haplotypes
        if(correction)
        {
            cout << "   Check: 3rd clustering... " << endl;
            maCluster.clear();
            paCluster.clear();
            unCluster.clear();
            pollens_updated.clear();       
            //
            pollens_corrected[0].poSeq = matGT_corrected;
            pollens_corrected[0].poAll = patGT_corrected;
        }
        if(correction && !cluster_pollens(pollens_corrected, 
                                          &maCluster,
                                          &paCluster,
                                          &unCluster,
                                          &pollens_updated)
          )              
        {
            return false;
        }
        // sub step 7. output phased markers along the current contig
        /*
            matGT_corrected again
            patGT_corrected again, being the final haplotypes
        */
        vector<unsigned long>::iterator mpositr     = mkrpos.begin();
        vector<unsigned long>::iterator mpositr_end = mkrpos.end();
        unsigned long iip = 0;
        while(mpositr != mpositr_end)
        {
            phasedmkrofp << "corr\t" 
                         << ctgid                          << "\t" << *mpositr                         << "\t"
                         << matGT_corrected.substr(iip, 1) << "\t" << patGT_corrected.substr(iip, 1)   << "\t"     
                         << matGT.substr(iip, 1)           << "\t" << patGT.substr(iip, 1)             << "\t"
                         << endl; 
            //
            iip     ++;
            mpositr ++;
        }
        // sub step 8. output cluster info of pollens
        //
        mergeofp << "\n#----Clustered----\n#" << endl;
        // find P/M/U info for this set of pollen
        map<int, string> pollenPMtmp;
        // 1. maternal cluster
        unsigned long matotalCoveredMarker = 0;
        unsigned long matotalBadMarker     = 0;      
        mergeofp << "\n" << matGT_corrected << "\t0_m_mkr corrected---- " << endl;
        ppitr     = maCluster.begin();
        ppitr_end = maCluster.end();
        while(ppitr != ppitr_end)
        {
            int pollenid = (*ppitr).first;
            POLLEN ptmp = (*ppitr).second;            
            mergeofp << ptmp.poSeq << "\t" << (*ppitr).first << "_x\t0_m" << endl;
            mergeofp << ptmp.poPat << "\t" << (*ppitr).first << "_x\t0_m" << endl;
            //
            matotalBadMarker     += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'P'); // should not be in M-cluster
            matotalCoveredMarker += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'P');
            matotalCoveredMarker += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'M');
            // PM=M
            //// string pmstring = get_crossover(ptmp.poPat, "M"); // check pollen in M-cluster
            
            
            
            unsigned long leftp = 0;
            double        score = 0;
            string     pmstring = get_break_pos(ptmp.poPat, "M", pollenid, &leftp, &score, ctgid, contigmarker);
            
            // output bed information if good
            bool outbed = false;
            vector<string> bedregions;
            if(!find_bed(ptmp.poPat, "M", pollenid, &outbed, &bedregions))
            {
                return false;
            }
            if(outbed)
            {
                vector<string>::iterator beditr;
                vector<string>::iterator beditr_end;
                beditr     = bedregions.begin();
                beditr_end = bedregions.end();
                while(beditr != beditr_end)
                {
                    vector<string> bedinfo    = split_string(*beditr, ':'); 
                    int leftindex             = atoi(bedinfo[1].c_str());
                    int rightindex            = atoi(bedinfo[2].c_str());
                    unsigned long leftmarker  = recover_real_marker_pos(leftindex,  contigmarker);    
                    unsigned long rightmarker = recover_real_marker_pos(rightindex, contigmarker); 
                    if(beditr == bedregions.begin()) 
                    {
                        leftmarker  = 1; // left-extending of the first region
                        if(ctgid.compare("tig00003688_pilonsB") == 0 )
                        {
                            leftmarker = 2131589;
                        }                        
                    }
                    if(rightindex == contigmarker.size() )
                    {
                        map<string, unsigned long>::iterator sizeitr = contigsize.find(ctgid);
                        if(sizeitr == contigsize.end())
                        {
                            cout << "   Erorr: " << ctgid << " not found in current size list. " << endl;
                            return false;
                        }
                        if(rightmarker > (*sizeitr).second)
                        {
                            cout << "   Error: markers propably derived from different contigs?" << endl;
                            return false;
                        }
                        rightmarker = (*sizeitr).second;
                    }                   
                    //                                                   
                    bedofp << ctgid       << "\t" 
                           << leftmarker  << "\t"
                           << rightmarker << "\t"
                           << bedinfo[1]  << "\t" 
                           << bedinfo[2]  << "\t" 
                           << bedinfo[0]  << "\t"
                           << "0_m"       << "\t"
                           << pollenid    << "_x" << endl;                    
                    beditr ++;
                }
            }
            
            pollenPMtmp.insert(std::pair<int, string>((*ppitr).first, pmstring)); // 2-bits: MP/PM/MM/PP
            ppitr ++;
        }
        // marker stat in M-cluster for current contig
        statmkrofp << ctgid                   << "\t"
                   << matGT_corrected.size()  << "\t"
                   << "M-cluster"             << "\t"
                   << maCluster.size()        << "\t"
                   << matotalCoveredMarker    << "\t"
                   << matotalBadMarker        << "\t"
                   << matotalCoveredMarker*1.0/matGT_corrected.size() << endl;
        // 2. paternal cluster
        unsigned long patotalCoveredMarker = 0;
        unsigned long patotalBadMarker     = 0;        
        mergeofp << "\n" << patGT_corrected << "\t0_p_mkr corrected---- " << endl;
        ppitr     = paCluster.begin();
        ppitr_end = paCluster.end();
        while(ppitr != ppitr_end)
        {
            int pollenid = (*ppitr).first;            
            POLLEN ptmp = (*ppitr).second;
            mergeofp << ptmp.poSeq << "\t" << (*ppitr).first << "_x\t0_p" << endl;
            mergeofp << ptmp.poPat << "\t" << (*ppitr).first << "_x\t0_p" << endl;
            //
            patotalBadMarker     += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'M'); // should not be in P-cluster
            patotalCoveredMarker += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'M');
            patotalCoveredMarker += std::count(ptmp.poPat.begin(), ptmp.poPat.end(), 'P');            
            // PM=P
            ////string pmstring = get_crossover(ptmp.poPat, "P"); // check pollen in P-cluster  
            
            
            
            unsigned long leftp = 0;
            double        score = 0;
            string pmstring     =  get_break_pos(ptmp.poPat, "P", pollenid, &leftp, &score, ctgid, contigmarker);  
            
               
            // output bed information if good
            bool outbed = false;
            vector<string> bedregions;
            if(!find_bed(ptmp.poPat, "P", pollenid, &outbed, &bedregions))
            {
                return false;
            }
            if(outbed)
            {
                vector<string>::iterator beditr;
                vector<string>::iterator beditr_end;
                beditr     = bedregions.begin();
                beditr_end = bedregions.end();
                while(beditr != beditr_end)
                {
                    vector<string> bedinfo    = split_string(*beditr, ':'); 
                    int leftindex             = atoi(bedinfo[1].c_str());
                    int rightindex            = atoi(bedinfo[2].c_str());
                    unsigned long leftmarker  = recover_real_marker_pos(leftindex,  contigmarker);    
                    unsigned long rightmarker = recover_real_marker_pos(rightindex, contigmarker);   
                    if(beditr == bedregions.begin()) 
                    {
                        leftmarker  = 1; // left-extending of the first region
                        if(ctgid.compare("tig00003688_pilonsB") == 0 )
                        {
                            leftmarker = 2131589;
                        }                        
                    }
                    if(rightindex == contigmarker.size() )
                    {
                        map<string, unsigned long>::iterator sizeitr = contigsize.find(ctgid);
                        if(rightmarker > (*sizeitr).second)
                        {
                            cout << "   Error: markers propably derived from different contigs?" << endl;
                            return false;
                        }
                        rightmarker = (*sizeitr).second;
                    } 
                    //                                                   
                    bedofp << ctgid       << "\t" 
                           << leftmarker  << "\t"
                           << rightmarker << "\t"
                           << bedinfo[1]  << "\t" 
                           << bedinfo[2]  << "\t" 
                           << bedinfo[0]  << "\t"
                           << "0_p"       << "\t"
                           << pollenid    << "_x" << endl;                    
                    beditr ++;
                }
            }               
               
                      
            pollenPMtmp.insert(std::pair<int, string>((*ppitr).first, pmstring)); // 2-bits: MP/PM/MM/PP
            ppitr ++;
        }
        // marker stat in M-cluster for current contig
        statmkrofp << ctgid                   << "\t"
                   << patGT_corrected.size()  << "\t"
                   << "P-cluster"             << "\t"
                   << paCluster.size()        << "\t"
                   << patotalCoveredMarker    << "\t"
                   << patotalBadMarker        << "\t"
                   << patotalCoveredMarker*1.0/patGT_corrected.size() << endl;
        // marker stat in P&M-cluster for current contig 
        unsigned long totalPollen  = maCluster.size()     + paCluster.size();
        unsigned long totalCovered = patotalCoveredMarker + matotalCoveredMarker;
        unsigned long totalBad     = patotalBadMarker     + matotalBadMarker;
        statmkrofp << ctgid                   << "\t"
                   << patGT_corrected.size()  << "\t"
                   << "P&M-cluster"           << "\t"
                   << totalPollen             << "\t"
                   << totalCovered            << "\t"
                   << totalBad                << "\t"
                   << totalCovered*1.0/patGT_corrected.size() << endl;                  
        // 3. un-clustered
        mergeofp << "\n#----un-clustered----\n" << endl;
        ppitr     = unCluster.begin();
        ppitr_end = unCluster.end();
        while(ppitr != ppitr_end)
        {
            int pollenid = (*ppitr).first;              
            POLLEN ptmp = (*ppitr).second;
            mergeofp << ptmp.poSeq << "\t" << (*ppitr).first << "_x\t0_u" << endl;
            mergeofp << ptmp.poPat << "\t" << (*ppitr).first << "_x\t0_u" << endl;
            // PM=U
//          pollenPMtmp.insert(std::pair<int, string>((*ppitr).first, "U" ));      
            // output bed information if good
            bool outbed = false;
            vector<string> bedregions;
            if(!find_bed(ptmp.poPat, "U", pollenid, &outbed, &bedregions))
            {
                return false;
            }
            if(outbed)
            {
                vector<string>::iterator beditr;
                vector<string>::iterator beditr_end;
                beditr     = bedregions.begin();
                beditr_end = bedregions.end();
                while(beditr != beditr_end)
                {
                    vector<string> bedinfo    = split_string(*beditr, ':'); 
                    int leftindex             = atoi(bedinfo[1].c_str());
                    int rightindex            = atoi(bedinfo[2].c_str());
                    unsigned long leftmarker  = recover_real_marker_pos(leftindex,  contigmarker);    
                    unsigned long rightmarker = recover_real_marker_pos(rightindex, contigmarker);   
                    if(beditr == bedregions.begin()) 
                    {
                        leftmarker  = 1; // left-extending of the first region
                        if(ctgid.compare("tig00003688_pilonsB") == 0 )
                        {
                            leftmarker = 2131589;
                        }                        
                    }
                    if(rightindex == contigmarker.size() )
                    {
                        map<string, unsigned long>::iterator sizeitr = contigsize.find(ctgid);
                        if(rightmarker > (*sizeitr).second)
                        {
                            cout << "   Error: markers propably derived from different contigs?" << endl;
                            return false;
                        }
                        rightmarker = (*sizeitr).second;
                    }                    
                    //                                                   
                    bedofp << ctgid       << "\t" 
                           << leftmarker  << "\t"
                           << rightmarker << "\t"
                           << bedinfo[1]  << "\t" 
                           << bedinfo[2]  << "\t" 
                           << bedinfo[0]  << "\t"
                           << "0_u"       << "\t"
                           << pollenid    << "_x" << endl;                    
                    beditr ++;
                }
            }                  
            pollenPMtmp.insert(std::pair<int, string>((*ppitr).first, "UU")); // 2-bits: UU            
            ppitr ++;
        }
        mergeofp << endl;  
        //
        mergeofp.close();
        //
        bedofp.close();
        // collect P/M/U info for this set of pollen
        contigsPM.insert(std::pair<string, map<int, string> >(ctgid, pollenPMtmp)); 
        //
        pcitr ++;
    }
    //
    // step 2: get contig genotype sequence at pollen samples; will be used to find 'linkage groups'
    map<string, BIPAT> contigPollenGTSeq;
    //
    gtinfo.str("");
    gtinfo << tmpfolder
           << "/s2_genotype_contig_seq.txt\0"; //
    ofstream ctgPMofp;
    ctgPMofp.open((gtinfo.str()).c_str(), ios::out);
    if(!ctgPMofp.good())
    {
        cout   << "   Error: cannot open file " << gtinfo.str() << endl;
        return false;
    }     
    map<string, map<int, string> >::iterator myctgitr;
    map<string, map<int, string> >::iterator myctgitr_end;
    myctgitr     = contigsPM.begin();
    myctgitr_end = contigsPM.end();
    while(myctgitr != myctgitr_end)
    {
        string ctgid                 = (*myctgitr).first;
        //
        string ctgPolGTSeq           = "";
        string leftGateSeq           = "";
        string rightGateSeq          = ""; 
        string recombined            = "";       
        map<int, string> pollenPMtmp = (*myctgitr).second;
        map<int, string>::iterator poitr;
        map<int, string>::iterator poitr_end;
        poitr     = pollenPMtmp.begin();
        poitr_end = pollenPMtmp.end();
        while(poitr != poitr_end)
        {
            //// ctgPMofp << (*poitr).second;    // 2-bits
            //// ctgPolGTSeq += (*poitr).second; // 2-bits
            string pmstring = (*poitr).second;
            assert(pmstring.size()==2);
            leftGateSeq    += pmstring.substr(0, 1);
            rightGateSeq   += pmstring.substr(1, 1);
            if(leftGateSeq.compare("u")  != 0 &&
               rightGateSeq.compare("u") != 0 &&
               pmstring.substr(0, 1).compare( pmstring.substr(1, 1) ) != 0)
            {
                recombined     += "R";
            }
            else
            {
                recombined     += "|";
            }
            poitr ++;
        }
        ////ctgPMofp << "\t";
        //
        ctgPMofp << leftGateSeq  << "\t" << ctgid << "\tleft-sequence"  << endl;  
        ctgPMofp << recombined   << "\t" << ctgid << "\tCO-info"        << endl;          
        ctgPMofp << rightGateSeq << "\t" << ctgid << "\tright-sequence" << endl;   
        // collect
        BIPAT tmpbipat;
        tmpbipat.leftPat  = leftGateSeq;
        tmpbipat.rightPat = rightGateSeq;
        contigPollenGTSeq.insert(std::pair<string, BIPAT>(ctgid, tmpbipat));     
        //
        myctgitr ++;
    }
    ctgPMofp.close();
    //
    // step 3: LG analysis
    if(!order_linked_contigs(contigPollenGTSeq, tmpfolder) )
    {
        return false;
    } 
    //   
    // s4 file where phased markers are collected
    phasedmkrofp.close();
    // s5 file where marker statistics along contigs are collected
    statmkrofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return 0;
}
//
unsigned long recover_real_marker_pos(int markerposition, map<unsigned long, ALLELE> contigmarker)
{
    // markerposition is 1-based positioning in pm pattern
    unsigned long realpos;
    //
    map<unsigned long, ALLELE>::iterator mitr;
    map<unsigned long, ALLELE>::iterator mitr_end;
    mitr     = contigmarker.begin();
    mitr_end = contigmarker.end();
    int i = 0;
    while(mitr != mitr_end)
    {   
        i ++;
        if(markerposition == i)
        {
            realpos = (*mitr).first;
        }
        mitr ++;
    }
    return realpos;
}     
// 
bool find_bed_v2(string poPat, string pmflag, int pollenid, bool* out, vector<string>* bedregions)
{
    /*
      PPuPPMMuMMuuuMMPuuuPuuMuM
      =>
      PPuPP:1:5
           MMuMMuuuMM:6:15
                     PuuuPuu:16:22
                            MuM:23:25
    */
    if(poPat.size()==0) return false; // should never happen
    //
    size_t posp = poPat.find("P");
    size_t posm = poPat.find("M");
    //
    if(posp==std::string::npos && posm==std::string::npos)
    {
        // only "U"
        *out = false;
    }else
    if(posp==std::string::npos && posm!=std::string::npos)
    {
        // only "M"
        std::stringstream ss;
        ss.str("");        
        ss << "M:" << "1:" << poPat.size(); // real position index now
        (*bedregions).push_back(ss.str());
        *out = true;
    }else
    if(posp!=std::string::npos && posm==std::string::npos)
    {
        // only "P"
        std::stringstream ss;
        ss.str("");               
        ss << "P:" << "1:" << poPat.size(); // real position index now
        (*bedregions).push_back(ss.str());      
        *out = true;
    }else
    {
        // "P&M&u"
        
        *out = true;
    }
    //   
    return true;
}                
//
bool find_bed(string poPat, string pmflag, int pollenid, bool* out, vector<string>* bedregions)
{
    /*
      PPuPPMMuMMuuuMMPuuuPuuMuM
      =>
      PPuPP:1:5
           MMuMMuuuMM:6:15
                     PuuuPuu:16:22
                            MuM:23:25
    */
    if(poPat.size()==0) return false; // should never happen
    //
    size_t posp = poPat.find("P");
    size_t posm = poPat.find("M");
    //
    if(posp==std::string::npos && posm==std::string::npos)
    {
        // only "U"
        *out = false;
    }else
    if(posp==std::string::npos && posm!=std::string::npos)
    {
        // only "M"
        std::stringstream ss;
        ss.str("");        
        ss << "M:" << "1:" << poPat.size(); // real position index now
        (*bedregions).push_back(ss.str());
        *out = true;
    }else
    if(posp!=std::string::npos && posm==std::string::npos)
    {
        // only "P"
        std::stringstream ss;
        ss.str("");               
        ss << "P:" << "1:" << poPat.size(); // real position index now
        (*bedregions).push_back(ss.str());      
        *out = true;
    }else
    {    
        // posp!=std::string::npos && posm!=std::string::npos
        // mixture of PMu
        string thisGT            = ""; // starting with P or M?
        size_t thisbegin         = 0;  // 0 based system - this needs to be increased by 1 to get real position index
        size_t thisend           = 0;  // 0 based system - this needs to be increased by 1 to get real position index
        size_t thisend_exclude_u = 0;
        if(posp < posm)
        {
            thisend            = posp;
            thisGT             = "P";
            thisend_exclude_u  = posp;                         
        }else
        {
            thisend            = posm;
            thisGT             = "M";
            thisend_exclude_u  = posm;            
        }
        for(int jj=thisend+1; jj<poPat.size(); jj ++)
        {
            if(poPat.substr(jj, 1).compare("u")==0) 
            {
                thisend           += 1; // extending 
                // reach end of pattern
                if(thisend==poPat.size()-1)
                {
                    // collect current region
                    std::stringstream ss;
                    ss.str("");  
                    ss << thisGT << ":" << thisbegin+1 << ":" << thisend+1; // real position index now        
                    (*bedregions).push_back(ss.str());
                    ss.str().clear();                       
                }
                continue;
            }else
            if(poPat.substr(jj, 1).compare(thisGT)==0)
            {
                thisend           += 1; // extending
                thisend_exclude_u  = jj;
                // reach end of pattern
                if(thisend==poPat.size()-1)
                {
                    // collect current region
                    std::stringstream ss;
                    ss.str("");  
                    ss << thisGT << ":" << thisbegin+1 << ":" << thisend+1; // real position index now        
                    (*bedregions).push_back(ss.str());
                    ss.str().clear();
                }
            }else
            {   // find a switch from thisGT to anotherGT (P->M or M->P)
                // collect current region
                std::stringstream ss;
                ss.str("");
                size_t left_marker_pos = floor((thisend_exclude_u+jj)/2);           // select the middle-left of the ambiguity region
                ss << thisGT << ":" << thisbegin+1 << ":" << left_marker_pos+1;     // real position index now
                (*bedregions).push_back(ss.str());
                ss.str().clear();
                // next region
                if(thisGT.compare("M") == 0) 
                {
                    thisGT = "P";
                }else
                if(thisGT.compare("P") == 0) 
                {
                    thisGT = "M";
                }else ;       
                /* 0-based marker sequence index, left_marker_pos+1 =< jj
                   theoretically, thisbegin is the middle marker position where PuuM or MuuuP happens. 
                      e.g, MuMuuuPuP
                           123456789
                      thisend_exclude_u = 3
                      jj                = 7
                      left_marker_pos = (3+7)/2 = 5 
                      next begin would be left_marker_pos+1 = 6 with new thisend_exclude_u = 6+1
                */
                thisbegin         = left_marker_pos+1;      
                thisend           = jj;
                thisend_exclude_u = jj;
            }
        }
        *out = true;
    }
    //
    return true;
}
// detect crossover along a pollen P/M pattern
string get_crossover(string poPat, string pmflag)
{
    /*
       not used after 20191101
       poPat  : 'P'/'M' pattern at markers along a contig: PPuPPMMuMMuuuMMMuuuuuMMuM
       pmflag : 'P' or 'M' to indicate cluster of current contig       
    */
    //
    // initialization
    string pmstring = "";
    //    
    string poPatbk = poPat;
    poPat.erase(std::remove(poPat.begin(), poPat.end(), 'u'), poPat.end()); // remove non-effective marekers with value: 'u'
    if(poPat.size()==0) 
    {
        pmstring = "UU";
        return pmstring;
    }
    //
    if(pmflag.compare("P") == 0)
    {
        // get count of 'M' in the beginning of poPat: to check whether it is MP-type pattern
        int Mcnt = 0;
        int Pcnt = 0;
        for(int li=0; li<poPat.size(); li++)
        {
            if(poPat.substr(li, 1).compare("M")==0)
            {
                Mcnt ++;
            }
            if(poPat.substr(li, 1).compare("P")==0)
            {
                Pcnt ++;
            }
            //            
            if(poPat.substr(li, 1).compare("P")==0 &&
               li >= poPat.size()/4) // caution on 1/4 cutoff
            {
                // only check the 1st 1/4 letters
                break;
            }
        }
        // first bit of pmstring
        if( poPat.substr(0, 1).compare("M")==0 && 
            Pcnt<Mcnt                          &&
           (Mcnt>=5 || ( poPat.size()<5 && Mcnt>0 && Mcnt>=poPat.size()/2) ) )// caution on 5   cutoff
        {
            pmstring += "M";
        }else
        {
            pmstring += "P";
        }
        // get count of 'M' in the end of poPat: to check whether it is PM-type pattern
        Mcnt = 0;
        Pcnt = 0;        
        for(int li=poPat.size()-1; li>=0; li--)
        {
            if(poPat.substr(li, 1).compare("M")==0)
            {
                Mcnt ++;
            }
            if(poPat.substr(li, 1).compare("P")==0)
            {
                Pcnt ++;
            } 
            if(poPat.substr(li, 1).compare("P")==0 &&
               li <= poPat.size()*3/4) // caution on 3/4 cutoff
            {
                // only check the last 1/4 letters
                break;
            }
        }
        // second bit of pmstring         
        if(poPat.substr(poPat.size()-1, 1).compare("M")==0 && 
           Pcnt<Mcnt                                       &&           
           ( Mcnt>=5 || (poPat.size()<5 && Mcnt>0 && Mcnt>=poPat.size()/2) ) ) // caution on 5   cutoff
        {
            pmstring += "M";
        }else
        {
            pmstring += "P";
        }     
        //
        cout << "   Check2rm: " << poPatbk << " ---- determined as " << pmstring << " for pollen in cluster P" << endl;
    }
    //
    if(pmflag.compare("M") == 0)
    {
        // get count of 'P' in the beginning of poPat: to check whether it is xM-type pattern
        int Pcnt = 0;
        int Mcnt = 0;
        for(int li=0; li<poPat.size(); li++)
        {
            if(poPat.substr(li, 1).compare("P")==0)
            {
                Pcnt ++;
            }
            if(poPat.substr(li, 1).compare("M")==0)
            {
                Mcnt ++;
            }
            if(poPat.substr(li, 1).compare("M")==0 &&
               li >= poPat.size()/4) // caution on 1/4 cutoff
            {
                // only check the 1st 1/4 letters
                break;
            }
        }
        // first bit of pmstring: if number of markers too small, the left-most must be "P" to determine it as "P"
        if(poPat.substr(0, 1).compare("P")==0 && 
           Mcnt<Pcnt                          &&                      
           ( Pcnt>=5 || ( poPat.size()<5 && Pcnt>0 && Pcnt>=poPat.size()/2) ) )// caution on 5   cutoff
        {
            pmstring += "P";
        }else
        {
            pmstring += "M";
        }
        // get count of 'P' in the end of poPat: to check whether it is xP-type pattern
        Pcnt = 0;
        Mcnt = 0;        
        for(int li=poPat.size()-1; li>=0; li--)
        {
            if(poPat.substr(li, 1).compare("P")==0)
            {
                Pcnt ++;
            }
            if(poPat.substr(li, 1).compare("M")==0)
            {
                Mcnt ++;
            }            
            if(poPat.substr(li, 1).compare("M")==0 &&
               li <= poPat.size()*3/4) // caution on 3/4 cutoff
            {
                // only check the last 1/4 letters
                break;
            }
        }  
        // second bit of pmstring: if number of markers too small, the right-most must be "P" to determine it as "P"
        if(poPat.substr(poPat.size()-1, 1).compare("P")==0 && 
           Mcnt<Pcnt                                       &&        
           ( Pcnt>=5 || (poPat.size()<5 && Pcnt>0 && Pcnt>=poPat.size()/2) ) )// caution on 5   cutoff
        {
            pmstring += "P";
        }else
        {
            pmstring += "M";
        }
        //
        cout << "   Check2rm: " << poPatbk << " ---- determined as " << pmstring << " for pollen in cluster M" << endl;
    }
    //
    return pmstring;
}
//  
string get_break_pos(string           poPat, 
                     string          pmflag,
                     int           pollenid,
                     unsigned long*   leftp, 
                     double*          score,
                     string                     contigname,
                     map<unsigned long, ALLELE> contigmarker)
{
    // used after 20191101
    // poPat: is the pollen PM pattern sequence
    // leftp : would be the i-th   marker position along the contig
    // rightp: would be the i+1-th marker position along the contig
    // score : would be score for the reported co
    // minimum number of markers hard-coded: at least 3
    //
    string pmstring = "";
    string poPatbk  = poPat;
    poPatbk.erase(std::remove(poPatbk.begin(), poPatbk.end(), 'u'), poPatbk.end()); // remove non-effective marekers with value: 'u'
    if(poPatbk.size()==0) 
    {
        pmstring = "UU";                // not sure on co
        *leftp   = 0;
        *score   = 0;
        cout << "   Check2rm: " << poPat << " ---- determined as "  << pmstring << " for pollen " << pollenid << "_x\t0_" << pmflag 
             << " (no effective markers) "                          << endl; 
        return pmstring;
    }
    if(poPatbk.size()==1)
    {
        if(poPatbk.compare("M")==0)
        {
            pmstring = "MM";
        }else
        if(poPatbk.compare("P")==0)
        {
            pmstring = "PP";
        }else ;
        *leftp   = 0; // no co
        *score   = 0; // no co     
        cout << "   Check2rm: " << poPat  << " ---- determined as "  << pmstring << " for pollen " << pollenid << "_x\t0_" << pmflag 
             << " ( 1 effective marker) "                            << endl; 
        return pmstring;  
    }
    //
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
            *leftp  = mi;  // (*leftp)-th marker is the left break marker, (*leftp+1) is the right break marker; 1-based positioning
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
    // 
    string checkfrom = "";
    // control on marker number: MuMMMMuuuuuuPuPPPu
    if(*score > 0 && poPatbk.size()<10)
    {
        // from cluster "M"; 
        if(pmflag.compare("M")    == 0)
        {
            // found "P" on left end;            
            if(pmstring.compare("PM") == 0)
            { 
               // at least 1 effective "P"; at least 1/3 "P"; first position was "P"
               if((poPatbk.substr(0, 1)).compare("P") == 0 && 
                  c1Lef > 0 && c1Lef >= poPatbk.size()/3.0)
                {
                    checkfrom = "1";
                }else
                {
                    pmstring  = "MM";
                    *leftp    = 0; // no co    
                    //*score  = 0; // no co   
                    checkfrom = "2";
                }
            }
            // found "P" on right end
            if(pmstring.compare("MP") == 0)
            {
               // at least 1 effective "P"; at least 1/3 "P"; last position was "P"                
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("P") == 0 && 
                   c1Rig > 0 && c1Rig >= poPatbk.size()/3.0)
                {
                    checkfrom = "3";
                }else
                {
                    pmstring = "MM";
                    *leftp   = 0; // no co
                    //*score  = 0; // no co 
                    checkfrom = "4";
                }
            }
        }
        // from cluster "P" 
        if(pmflag.compare("P")    == 0)
        {
            // but found "M" on left end
            if(pmstring.compare("MP") == 0 )
            {            
               // at least 1 effective "M"; at least 1/3 "M"; first position was "M"                
                if((poPatbk.substr(0, 1)).compare("M") == 0 && 
                   c2Lef > 0 && c2Lef >= poPatbk.size()/3.0)
                {
                    checkfrom = "5";
                }else
                {
                    pmstring = "PP";
                    *leftp   = 0; // no co
                    //*score  = 0; // no co    
                    checkfrom = "6";
                }
            }
            // but found "M" on right end
            if(pmstring.compare("PM") == 0)
            {            
               // at least 1 effective "M"; at least 1/3 "M"; last position was "M"                                
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("M") == 0 &&     
                   c2Rig > 0 && c2Rig >= poPatbk.size()/3.0)
                {
                    checkfrom = "7";
                }else
                {
                    pmstring = "PP";
                    *leftp   = 0; // no co
                    //*score  = 0; // no co 
                    checkfrom = "8";
                }
            }
        }                
        //
    }else
    if(*score>=minCOscore && poPatbk.size()>10)
    {
        // from cluster "M" 
        if(pmflag.compare("M")    == 0)
        {        
            // but found "P" on left end 
            if(pmstring.compare("PM") == 0)
            {            
                // and at least 5 effective "P" and first position was "P"
                if((poPatbk.substr(0, 1)).compare("P") == 0 && c1Lef >= 5 )
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
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("P") == 0 && c1Rig >= 5 )
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
        if(pmflag.compare("P")    == 0)
        {        
            // but found "M" on left end 
            if(pmstring.compare("MP") == 0)
            {            
                // and at least 5 effective "M" and first position was "M"                
                if((poPatbk.substr(0, 1)).compare("M") == 0 && c2Lef >= 5)
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
            // from cluster "P" but found "M" on left end and effective last  position was "M"
            if(pmstring.compare("PM") == 0)
            {            
                if((poPatbk.substr(poPatbk.size()-1, 1)).compare("M") == 0 && c2Rig >= 5)
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
        if(pmflag.compare("M") == 0)
        {
            pmstring     = "MM";
            checkfrom = "16";
        }else
        if(pmflag.compare("P") == 0)
        {
            pmstring     = "PP";
            checkfrom = "17";
        }
        *leftp   = 0; // no co
        //*score = 0; // no co  
    }
    //
    cout << "   Check2rm: " << poPat  << " ---- determined as "  << pmstring << " for pollen " << pollenid << "_x\t0_" << pmflag 
         << " (CO after "   << *leftp << "th marker: score "     << *score   << " - case "     << checkfrom << " ) "; // *** 
    // output co interval
    if(*leftp != 0)
    {
        int leftindex  = *leftp-1;
        while(leftindex>0)
        {
            if(poPat.substr(leftindex, 1).compare("u") == 0)
            {
                leftindex --;
            }else
            {
                break;
            }
        }
        int rightindex = *leftp;
        while(rightindex<poPat.size()-1)
        {
            if(poPat.substr(rightindex, 1).compare("u") == 0)
            {
                rightindex ++;
            }else
            {
                break;
            }
        }
        // output CO real positions of markers: increase to 1-based positioning <= leftindex+1; real position in poPat
        unsigned long leftmarker  = recover_real_marker_pos(leftindex+1,  contigmarker);
        unsigned long rightmarker = recover_real_marker_pos(rightindex+1, contigmarker);
        cout << " (CO interval " 
             << contigname 
             << " "
             << leftindex+1  << "th.mkr=" << leftmarker 
             << "-"
             << rightindex+1 << "th.mkr=" << rightmarker
             << ", size="
             << rightmarker - leftmarker + 1 << ")"
             << endl;
    }else
    {
        cout <<  endl; // to end ***
    }
    return pmstring;
}
//
bool correct_haplotype3_with_pmClusterCnt(map<int, POLLEN> maCluster,
                                          map<int, POLLEN> paCluster,
                                          string*          matGT,
                                          string*          patGT)
{
    // this is stage 3 correction: use 'P/M' counts in Paternal and Maternal clusters to correct point phasing    
    /*
        1) WITHOUT CORRECTION:    
             Paternal cluster at position i: P=50, M=02, u=200
             Maternal cluster at position i: P=03, M=45, u=250           
        2) WITH    CORRECTION:
             Paternal cluster at position i: P=02, M=50, u=200
             Maternal cluster at position i: P=45, M=03, u=250  
        
        struct CLUSCNT
        {
            int Pcnt;
            int Mcnt;
            int ucnt;
        };               
    */
    if((*patGT).size() != (*matGT).size() )
    {
        cout << "   Error: inconsistency in haplotype lengths. " << endl;
        return false;
    } 
    unsigned long markerNum = (*patGT).size();
    // 0. initialization on counts for two clusters
    map<int, CLUSCNT> maClusterCnt;
    map<int, CLUSCNT> paClusterCnt;        
    for(int ci = 0; ci < markerNum; ci ++)
    {
        CLUSCNT posCNT;
        posCNT.Pcnt = 0;
        posCNT.Mcnt = 0;
        posCNT.ucnt = 0;
        paClusterCnt.insert(std::pair<int, CLUSCNT>(ci, posCNT));
        maClusterCnt.insert(std::pair<int, CLUSCNT>(ci, posCNT));        
    }
    //
    map<int, POLLEN>::iterator ppitr;
    map<int, POLLEN>::iterator ppitr_end;
    map<int, CLUSCNT>::iterator cntitr;
    // 1. check maternal cluster attached to 0_m_mkr: matGT
    ppitr     = maCluster.begin();
    ppitr_end = maCluster.end();
    while(ppitr != ppitr_end)
    {
        POLLEN ptmp         = (*ppitr).second;
        string this_pattern = ptmp.poPat; // MMMMPPMMMMMMuM...MMMuMMMM
        // check P/M/u
        for(int ci = 0; ci < markerNum; ci ++)
        {
            // get "P" or "M" at current position
            string posPattern = this_pattern.substr(ci, 1);
            CLUSCNT posCNT;
            get_PM_count(posPattern, &posCNT);            
            // 
            cntitr = maClusterCnt.find(ci);
            // update
            (*cntitr).second.Pcnt += posCNT.Pcnt;
            (*cntitr).second.Mcnt += posCNT.Mcnt;
            (*cntitr).second.ucnt += posCNT.ucnt;
        }
        //
        ppitr ++;
    }
    // 2. check paternal cluster attached to 0_p_mkr: patGT
    ppitr     = paCluster.begin();
    ppitr_end = paCluster.end();
    while(ppitr != ppitr_end)
    {
        POLLEN ptmp         = (*ppitr).second;
        string this_pattern = ptmp.poPat; // PPPPMMPPPuPPPP...PPPPPuPP
        // check P/M/u
        for(int ci = 0; ci < markerNum; ci ++)
        {
            // get "P" or "M" at current position
            string posPattern = this_pattern.substr(ci, 1);
            CLUSCNT posCNT;
            get_PM_count(posPattern, &posCNT);            
            // 
            cntitr = paClusterCnt.find(ci);
            // update
            (*cntitr).second.Pcnt += posCNT.Pcnt;
            (*cntitr).second.Mcnt += posCNT.Mcnt;
            (*cntitr).second.ucnt += posCNT.ucnt;
        }
        //
        ppitr ++;
    }
    // 3. correct haplotypes according to counts
    for(int ci=0; ci<markerNum; ci ++)
    {
        // get counts on P/M in paternal cluster for current position
        map<int, CLUSCNT>::iterator pacitr = paClusterCnt.find(ci);
        assert(pacitr != paClusterCnt.end() );
        CLUSCNT paPosCNT = (*pacitr).second;
        // get counts on P/M in maternal cluster for current position
        map<int, CLUSCNT>::iterator macitr = maClusterCnt.find(ci);
        assert(macitr != maClusterCnt.end() );
        CLUSCNT maPosCNT = (*macitr).second;     
        // swap paternal/maternal haplotypes
        if( (paPosCNT.Pcnt<paPosCNT.Mcnt) && 
            (maPosCNT.Pcnt>maPosCNT.Mcnt)    )
        {
            string mattype = (*matGT).substr(ci, 1); // "P/M" at maternal pattern
            string pattype = (*patGT).substr(ci, 1); // "P/M" at paternal pattern
            (*matGT).replace(ci, 1, pattype);
            (*patGT).replace(ci, 1, mattype);
        }
    }  
    cout << "   Check: 3rd corrected haplotypes: " << endl
         << "\t"  << *matGT                        << endl
         << "\t"  << *patGT                        << endl;                  
    //    
    return true;
}                        
//
bool get_PM_count(string posPattern, CLUSCNT* posCNT)
{
    if(posPattern.compare("P") == 0)
    {
        (*posCNT).Pcnt = 1;
        (*posCNT).Mcnt = 0;
        (*posCNT).ucnt = 0;                          
    }else
    if(posPattern.compare("M") == 0)
    {
        (*posCNT).Pcnt = 0;
        (*posCNT).Mcnt = 1;
        (*posCNT).ucnt = 0;                          
    }else   
    {
        (*posCNT).Pcnt = 0;
        (*posCNT).Mcnt = 0;
        (*posCNT).ucnt = 1;                          
    }
    return true;
}

//
bool correct_haplotype2_with_consensus(
                         map<int, POLLEN>                   maCluster,
                         map<int, POLLEN>                   paCluster,
                         string*                            matGT,
                         string*                            patGT)
{
    // this is stage 2 correction: using consensus of all pollen to determine maternal haplotype - used
    // from flip_pcluster function, patenral pollen seqs have been converted to maternal seqs
    map<int, POLLEN> bothCluster = maCluster;
    bothCluster.insert(paCluster.begin(), paCluster.end());
    int pollenNum                = bothCluster.size();    
    
    // position-wise consensus
    for(unsigned long ii = 0; ii < (*matGT).size(); ii ++)
    {
        map<string, int> baseCnt;             // parental allele cnt
        map<string, int>::iterator baseitr;
        //
        map<int, POLLEN>::iterator ppitr;     // traverse each collected pollen
        map<int, POLLEN>::iterator ppitr_end; // traverse each collected pollen   
        ppitr     = bothCluster.begin();
        ppitr_end = bothCluster.end();
        while(ppitr != ppitr_end)
        {
            POLLEN ptmp  = (*ppitr).second;
            string poSeq = ptmp.poSeq;   
            
            string base  = poSeq.substr(ii, 1);// can be ACGT/N
            
            baseitr = baseCnt.find(base);
            if(baseitr != baseCnt.end())
            {
                (*baseitr).second ++;
            }else
            {
                baseCnt.insert(std::pair<string, int>(base, 1));
            }
            ppitr ++;
        }
        // find major allele and minor allele
        baseitr = baseCnt.begin();
        map<string, int>::iterator baseitr_end = baseCnt.end();
        string majBase = "";
        int    majCnt  = 0;
        string minBase = "";
        int    minCnt  = 0;
        while(baseitr != baseitr_end)
        {
            string base = (*baseitr).first;
            int    cnt  = (*baseitr).second;
            if(base.compare("A")==0 ||
               base.compare("C")==0 || 
               base.compare("G")==0 ||
               base.compare("T")==0  )
            {
                if(cnt >  majCnt) // to consider 0 count
                {
                    majBase = base;
                    majCnt  = cnt;
                }else
                if(cnt >= minCnt) // to consider 0 count
                {
                    minBase = base;
                    minCnt  = cnt;
                }
                else ;              
            }
            baseitr ++;
        }
        //
        string matBase = (*matGT).substr(ii, 1);  // the base at current maternal haplotype, which should differ paternal one
        if(majBase.compare(minBase) != 0)               
        if(matBase.compare(majBase) != 0)         // swap as main base in maternal
        {
            if(majCnt > 0)
            {
                (*matGT).replace(ii, 1, majBase); //
                (*patGT).replace(ii, 1, matBase); //
            }
        }
        // if equal base then no change in the orignal pattern
    }
    cout << "   Check: 2nd corrected haplotypes: " << endl
         << "\t"  << *matGT                        << endl
         << "\t"  << *patGT                        << endl;     
    return true;
} 
//
// flip pcluster to be consistent and merged with mcluster and find the consensus to get the new maternal haplotype
bool flip_pcluster(map<int, POLLEN>* paCluster, string matGT, string patGT)
{
    /*
      struct POLLEN
      {
          string poSeq;   // seq at markers for pollen haplotype
          string poAll;   // over all, can be determined as 'M', 'P', 'MP', 'PM', or 'N'
          string poPat;   // string like MMMMMMppMMMMMMMpMMMMMMMM: indicating parental genotypes
          vector<unsigned long> mkrpos; // position of markers for this contig
      };    
    */
    if((*paCluster).size()==0) return true;
    //
    map<int, POLLEN>::iterator ppitr;
    map<int, POLLEN>::iterator ppitr_end;
    ppitr     = (*paCluster).begin();
    ppitr_end = (*paCluster).end();
    while(ppitr != ppitr_end)
    {
        POLLEN ptmp = (*ppitr).second;
        
        string poSeq = ptmp.poSeq;
        string poAll = ptmp.poAll;
        string poPat = ptmp.poPat;
        
        for(unsigned long ii = 0; ii < poPat.size(); ii ++)
        {
            if(poPat.substr(ii, 1).compare("P") == 0)
            {
                // flipping pattern
                (*ppitr).second.poPat = "M";
                (*ppitr).second.poSeq.replace(ii, 1, matGT.substr(ii, 1));
            }else
            if(poPat.substr(ii, 1).compare("M") == 0)
            {
                // flipping pattern
                (*ppitr).second.poPat = "P";
                (*ppitr).second.poSeq.replace(ii, 1, patGT.substr(ii, 1));                
            }else ;           
        }
        //
        (*ppitr).second.poAll = "M";
                
        // next paternal pollen
        ppitr ++;
    }
    return true;
}
//
// correct raw haplotypes from assembly with pollen haplotypes                         
bool correct_haplotype(  map<int, POLLEN>                   pollens,
                         map<int, POLLEN>*                  pollens_corrected)
{
    // this is stage 1 correction: using pollen individuals to correct assembled m/paternal haplotypes -used
    /*
       for one contig:
       
       // parents
       AACTGATCTGTGTTCCAATAGCACACGGGGAGCACGGATTTGTTCGGACGTATCGAATGCA	0_m::matGT initial backbone
       GTGCATCGCACACCTTGGCGATGGTTCAAATATGTTAGCGATACTAAGTACGCTCTCGCAC	0_p::patGT initial backbone
       // pollens:
       NNNNNNNNNNNNNNNNNNNNANNNNNNNNNNNNNNTNNNNNNNNNNNNNNNNNNNNNGCCC	1_x
       NTGNNNNNNNCNNNNNNNCGANGGNNCNAATATGTTAGCGATNCTAAGTACGCTCTCGCAC	2_x
       NTGCATCNNNCNNNNNNNNNNNNNNNNNNNNNNNNTNNNNNNNNNNNNNNNNNNNNNNNCC	3_x
       GTGCATCGCACUCCTTGGCGATGGTTCANATNNGTTNNNNNNNNNNNNNNNNNNNNNNNNN	4_x
       NNNNNNNNNNNNNTCNAANNNCNNNNNNNNNNNANGNNNNNNNNNNNNNNNNNNNNNTGCA	5_x
       UUUCUUUUUUUUUCUUUUUUUUGUTTCAUAUUUGTUUUUUATACTAAGTACGCUUTCGCAC	6_x
       GUUCATCUUUCUUUUTGUUUUUUUUTUUUUUUUUUUUUUUUUUUUUUUUUUUUUGUUUUUU	7_x
       NNNNNNNNNNCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGNNN	8_x
       UANTGNNCTNNNNNNCNANNGNNNNNNAANNNNNCGGATTUGTNCGGACGNNTCGAAGUUU	9_x
       ANNNGATCTGUGTNNNNNTAGNACACNGNGNGCANGGATNTNNNCGNNNNNNNNNNNNNNN	10_x
       GTGCATCGCACUCCTTGGCGAUGGTTCAAATATGTTAGCGATACTAAGTACGCTCTCGCAC	11_x
       GUGCATCGCACUCCTTGGCGATGGTTCAAATATGTTAGCGATACTAAGTACGCTCTCGCUC	12_x
       GTGCATCGCACACCTTGGUUUTGGTTCUAAUUTGUTAGCGATACTAAGTACGCTCTCGCUC	13_x
       GNGCATCGCANACCTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUC	14_x
       NNNNNNNNNGTGTTCCAATAGCACACNGGGAGCANGGATTTGNNNNNNNNNNNNNNNNNNN	15_x    
    */
    // this is to correct maternal and paternal haplotypes, or phasing
    string matGT = pollens[0].poSeq; // initial maternal backbone
    string patGT = pollens[0].poAll; // initial paternal backbone
    // update result
    if(matGT.size() == 1 && patGT.size() == 1)
    {
        (*pollens_corrected)          = pollens;    
        (*pollens_corrected)[0].poSeq = matGT; // 1 bp maternal backbone requiring no correction
        (*pollens_corrected)[0].poAll = patGT; // 1 bp paternal backbone requiring no correction  
        return true;    
    }
    //
    string matGT_corrected = "";
    string patGT_corrected = "";
    // 
    // initialize dimers at each pair of adjacent snps
    map<unsigned long, map<string, int> > dimers;    // <position, <2snp-dimer, count> >
    //
    for(unsigned long ii = 0; ii < matGT.size() - 1; ii ++)
    {
        // get     4 alleles from 2 snps: e.g., snp 1: (A, G) and snp 2: (A, T)
        string pos1allele1 = matGT.substr(ii+0, 1);  // A
        string pos1allele2 = patGT.substr(ii+0, 1);  // G
        string pos2allele1 = matGT.substr(ii+1, 1);  // A
        string pos2allele2 = patGT.substr(ii+1, 1);  // T
        // combine 4 alleles from 2 snps
        map<string, int> dimers_at_pos;                
        string dimertmp = "";
        // dimer 1
        dimertmp = pos1allele1 + pos2allele1;        // AA
        dimers_at_pos.insert(std::pair<string, int>(dimertmp, 0));        
        // dimer 2
        dimertmp = pos1allele1 + pos2allele2;        // AT
        dimers_at_pos.insert(std::pair<string, int>(dimertmp, 0));             
        // dimer 3
        dimertmp = pos1allele2 + pos2allele1;        // GA
        dimers_at_pos.insert(std::pair<string, int>(dimertmp, 0));             
        // dimer 4
        dimertmp = pos1allele2 + pos2allele2;        // GT
        dimers_at_pos.insert(std::pair<string, int>(dimertmp, 0));             
        // collect 4 dimers at postion ii 
        dimers.insert(std::pair<unsigned long, map<string, int> >(ii, dimers_at_pos));
    }
    //
    map<int, POLLEN>::iterator ppitr;
    map<int, POLLEN>::iterator ppitr_end;
    ppitr     = pollens.begin();
    ppitr_end = pollens.end();
    while(ppitr != ppitr_end)
    {
        POLLEN ptmp     = (*ppitr).second;
        // get current pollen sequence at markers along contigs
        string pollenGT = ptmp.poSeq;
        // check consistency in sequence length
        if (matGT.size() != patGT.size()    || 
            matGT.size() != pollenGT.size() ||
            patGT.size() != pollenGT.size() ) 
        {
            cout << "   Error: sequences of pollen and haplotype initials are with different lengths. " << endl;
            return false;
        }
        // vote each initialized dimmer according parental snps
        map<unsigned long, map<string, int> >::iterator positr; 
        for(unsigned long jj = 0; jj < pollenGT.size() - 1; jj ++)
        {
            // get a dimer along current pollen sequence
            string dimertmp = pollenGT.substr(jj, 2);
            // find  initialized dimmers according to parental un-phased snps at the current position
            positr = dimers.find(jj);
            assert(positr != dimers.end());            
            // score initialized dimmers according to parental un-phased snps at the current position
            map<string, int>::iterator dimeritr;
            dimeritr = ((*positr).second).find(dimertmp);
            if(dimeritr != ((*positr).second).end())
            {
                (*dimeritr).second += 1;
            }
        }
        // next pollen
        ppitr ++;
    }
    // get corrected haplotypes
    map<unsigned long, map<string, int> >::iterator positr;
    map<unsigned long, map<string, int> >::iterator positr_end;
    positr         = dimers.begin();
    positr_end     = dimers.end();
    string lastMat = "";
    string lastPat = "";
    while(positr != positr_end)
    {
        map<string, int> dimers_at_pos;
        dimers_at_pos = (*positr).second;
        //
        map<string, int>::iterator ditr;
        map<string, int>::iterator ditr_end;
        // find mostly voted genotype
        ditr     = dimers_at_pos.begin();
        ditr_end = dimers_at_pos.end();
        string maxDimer     = "";
        int    maxCnt       = 0;
        string left_allele  = "";
        string right_allele = "";
        while(ditr != ditr_end)
        {
            string this_dimer = (*ditr).first; 
            //           
            string this_left  = this_dimer.substr(0, 1);
            if(left_allele.find(this_left) == std::string::npos)
            {
                left_allele += this_left;
            }
            string this_right  = this_dimer.substr(1, 1);
            if(right_allele.find(this_right) == std::string::npos)
            {
                right_allele += this_right;
            }
            //            
            if( (*ditr).second > maxCnt)
            {
                maxDimer = (*ditr).first;
                maxCnt   = (*ditr).second;
            }
            ditr ++;
        }
        if(maxCnt == 0) // not voted cases: there might be a few; randomly set the first as best
        {
            maxDimer = (*(dimers_at_pos.begin() ) ).first;
            maxCnt   = (*(dimers_at_pos.begin() ) ).second;            
        }
        // find complemented genotype to the genotype with maximum vote
        string max_left  = maxDimer.substr(0, 1);        
        string max_right = maxDimer.substr(1, 1);
        string cmp_left  = ""; // complemented left  allele
        string cmp_right = ""; // complemented right allele
        for (int kk = 0; kk < 2; kk ++)
        {
            if(left_allele.substr(kk, 1).compare(max_left) != 0) 
            {
                cmp_left = left_allele.substr(kk, 1);
            }
            if(right_allele.substr(kk, 1).compare(max_right) != 0) 
            {
                cmp_right = right_allele.substr(kk, 1);
            }
        }
        string cmpDimer = cmp_left + cmp_right;
        assert( dimers_at_pos.find(cmpDimer) != dimers_at_pos.end() );
        int    cmpCnt   = dimers_at_pos[cmpDimer];
        //
        // check the other possibility: total support for both combination of alleles.
        string dimer3 = cmpDimer.substr(0, 1) + maxDimer.substr(1, 1);
        string dimer4 = maxDimer.substr(0, 1) + cmpDimer.substr(1, 1);  
        assert(dimers_at_pos.find(dimer3) != dimers_at_pos.end());
        assert(dimers_at_pos.find(dimer4) != dimers_at_pos.end());
        int    cnt3   = dimers_at_pos[dimer3];
        int    cnt4   = dimers_at_pos[dimer4];
        bool mcorrected=false;
        if(cnt3+cnt4  > cmpCnt + maxCnt)
        {
            cmpCnt   = cnt3;
            cmpDimer = dimer3;
            maxCnt   = cnt4;
            maxDimer = dimer4;
            if(cmpCnt > maxCnt)
            {
                int tmpcnt    = cmpCnt;
                cmpCnt        = maxCnt;
                maxCnt        = tmpcnt;
                //
                string tmpstr = cmpDimer;
                cmpDimer      = maxDimer;
                maxDimer      = tmpstr;
                //
                mcorrected    = true;
            }
        }        
        //        
        // out mostly voted genotype
        if(lastMat.size () == 0)
        {
            cout << "\t" << maxDimer << ":" << maxCnt << "\t"
                 <<         cmpDimer << ":" << cmpCnt << "\tother\t";
            // collect phased haplotypes                 
            matGT_corrected = maxDimer;
            patGT_corrected = cmpDimer;  
            // update
           lastMat = maxDimer;
           lastPat = cmpDimer;                             
        }else
        {
            string lastMatEnd = lastMat.substr(1, 1); // end of last dimmer for maternal haplotype
            string lastPatEnd = lastPat.substr(1, 1); // end of last dimmer for paternal haplotype
            string thisSta1   = maxDimer.substr(0, 1);// sta of this dimmer with maximum vote 
            string thisSta2   = cmpDimer.substr(0, 1);// sta of this dimmer complemented maximum vote
            //
            if(lastMatEnd.compare(thisSta1) == 0 && 
               lastPatEnd.compare(thisSta2) == 0)
            {
                cout << "\t" << maxDimer << ":" << maxCnt << "\t"
                     <<         cmpDimer << ":" << cmpCnt << "\tother\t";  
                // collect phased haplotypes                 
                matGT_corrected += maxDimer.substr(1,1);
                patGT_corrected += cmpDimer.substr(1,1);    
                // update
                lastMat = maxDimer;
                lastPat = cmpDimer;                                                        
            }else
            {
                cout << "\t" << cmpDimer << ":" << cmpCnt << "\t"
                     <<         maxDimer << ":" << maxCnt << "\tother\t";    
                // collect phased haplotypes                 
                matGT_corrected += cmpDimer.substr(1,1);
                patGT_corrected += maxDimer.substr(1,1); 
                // update
                lastMat = cmpDimer;
                lastPat = maxDimer;                                                            
            } 
                              
        } 
        // check 20191025
        map<string, int>::iterator oditr;
        map<string, int>::iterator oditr_end;    
        oditr     = dimers_at_pos.begin();
        oditr_end = dimers_at_pos.end();       
        while(oditr != oditr_end)
        {
            string this_dimer = (*oditr).first; 
            if(this_dimer.compare(maxDimer)!=0 && this_dimer.compare(cmpDimer)!=0)
            {
                cout << this_dimer << ":" << (*oditr).second << "\t";
            }
            oditr ++;
        }               
        if(mcorrected)
        {
            cout << "\tcorrected";
        }  
        cout << endl;               
        // next position
        positr ++;
    }
    //
    cout << "   Check: initial   haplotypes: " << endl
         << "\t"  << matGT                     << endl
         << "\t"  << patGT                     << endl;
    cout << "   Check: corrected haplotypes: " << endl
         << "\t"  << matGT_corrected           << endl
         << "\t"  << patGT_corrected           << endl;     
    // update result
    (*pollens_corrected) = pollens;    
    (*pollens_corrected)[0].poSeq = matGT_corrected; // corrected maternal backbone
    (*pollens_corrected)[0].poAll = patGT_corrected; // corrected paternal backbone    
    //
    return true;
}

//
double calculate_lod(string patA, string patB)
{
    /* contigPollenGTSeq:
       PPPPMPPPMMPPPPM patA
       PPPMMPPMPMPPUMM patB
    */
    double lod = 0.0;
    assert(patA.size() == patB.size());
    double r   = 0.0; // count of     recombined case 
    double nr  = 0.0; // count of non-recombined case
    for(int ii=0; ii<patA.size(); ii++)
    {
        string Ageno = patA.substr(ii, 1);
        string Bgeno = patB.substr(ii, 1);
        if(Ageno.compare("U")  !=0 &&
           Bgeno.compare("U")  !=0 &&
           Ageno.compare(Bgeno)==0   )
        {
            nr ++;
        }else
        {
            r  ++;
        }
    }
    double totalEffective = r + nr;
    if(totalEffective > 0)
    {
        // recombination rate
        double theta          = r / totalEffective;
        //
        double thisNumerator;
        double thisDenominator;
        if(theta!=0  && theta!=1)
        {
            thisNumerator     = pow(theta, r) * pow(1-theta, nr);
        }else
        {
            thisNumerator     = 1;
        }
        thisDenominator       = pow(0.5, totalEffective);    
        //
        lod                   = log10(thisNumerator / thisDenominator);
    }else
    {
        lod                   = -1000.0; // cannot be determined
    }
    //
    return lod;
}
//
bool order_linked_contigs(map<string, BIPAT> contigPollenGTSeq, string tmpfolder)
{
    /*  contigPollenGTSeq:
	   PMPPMPPPMMPPPPM tig00000013_pilon left  sequence
	   PPPPMPPPMMPPPPM tig00000013_pilon right sequence
	   //
   	   PPPMMPPMPMPPUMM tig00000026_pilon ...
   	   PPPMMPPMPMPMUMM tig00000026_pilon ...  	  
   	   // 
	......
    */
    map<string, int> collected_edge; // if a->b edge exists, then ignore b->a edge
    // prepare dot output file
    string dotinfo = tmpfolder + "/s3_genotype_contig_GT_ordered_LGlike.dot\0"; // 
    ofstream dotofp;
    dotofp.open(dotinfo.c_str(), ios::out);
    if(!dotofp.good())
    {
        cout   << "   Error: cannot open file " << dotinfo << endl;
        return false;
    } 
    dotofp << "/* this graph is translated from similarity of genotype sequences of contigs */" << endl;
    dotofp << "digraph\tGraph_1 {" << endl;
    //
    map<string, BIPAT>::iterator ctgitr;
    map<string, BIPAT>::iterator ctgitr_end;
    ctgitr     = contigPollenGTSeq.begin();
    ctgitr_end = contigPollenGTSeq.end();
    while(ctgitr != ctgitr_end)
    {
        string this_contig_ID = (*ctgitr).first;
        string this_contig_GT = (*ctgitr).second.rightPat; // out vertex
        //
        CTGLG contigBestLods;
        contigBestLods.lod[0] = -1000; // best
        contigBestLods.lod[1] = -1000; // second-best
        contigBestLods.edge1  = "";
        contigBestLods.edge2  = "";
        // second level traverse of contigPollenGTSeq;
        map<string, BIPAT>::iterator ctgitr2;
        map<string, BIPAT>::iterator ctgitr2_end;
        ctgitr2     = contigPollenGTSeq.begin();
        ctgitr2_end = contigPollenGTSeq.end();
        while(ctgitr2 != ctgitr2_end)
        {
            string tmp2_contig_ID = (*ctgitr2).first;
            string tmp2_contig_GT = (*ctgitr2).second.leftPat; // in vertex
            int effectiveMkrNum   = 0;
            effectiveMkrNum      += std::count(tmp2_contig_GT.begin(), tmp2_contig_GT.end(), 'P');
            effectiveMkrNum      += std::count(tmp2_contig_GT.begin(), tmp2_contig_GT.end(), 'M');            
            // IF NEXT CONTIG DOES NOT HAVE ENOUGH PATERNAL/MATERNAL INFORMATION, SKIP IT 
            if(effectiveMkrNum < tmp2_contig_GT.size()/3)
            {
                 ctgitr2 ++;
                 continue;                 
            }
            //
            if(tmp2_contig_ID.compare(this_contig_ID) != 0)
            {
                // reversed edge not in collected list
                string edge1 = this_contig_ID + "->" + tmp2_contig_ID;
                string edge2 = this_contig_ID + "->" + tmp2_contig_ID + "_rev";
                if(collected_edge.find(edge1) == collected_edge.end() &&
                   collected_edge.find(edge2) == collected_edge.end() )
                {
                    /*
                    // original
                    int score1 = 0;
                    if(!get_match_score(this_contig_GT,        tmp2_contig_GT,       &score1) ) return false;
                    // reversed
                    string tmp2_contig_GT_reversed = "";
                    if(!reverse_contig_GT(tmp2_contig_GT, &tmp2_contig_GT_reversed) )           return false;
                    int score2 = 0;
                    if(!get_match_score(this_contig_GT,    tmp2_contig_GT_reversed,  &score2) ) return false;
                    // output as dot
                    if(score1 >= score2)
                    {
                        if(score1 >= this_contig_GT.size()*0.90 )
                        {
                            // select original: 	100052 -> 358538 [color=aquamarine, penwidth=1, arrowsize=1, label=45];
                            dotofp << "\t"   << this_contig_ID 
                                   << " -> " << tmp2_contig_ID
                                   << " "    << "[color=black, penwidth=1, arrowsize=1, label=" << score1 << "];"    
                                   << endl;
                            // collect edge
                            collected_edge.insert(std::pair<string, int>(this_contig_ID + "->" + tmp2_contig_ID, 1));       
                        }                    
                    }else
                    {
                        if(score2 >= this_contig_GT.size()*0.90 )
                        {                
                            // select reversed    
                            dotofp << "\t"   << this_contig_ID 
                                   << " -> " << tmp2_contig_ID << "_rev"
                                   << " "    << "[color=red, penwidth=1, arrowsize=1, label=" << score2 << "];"
                                   << endl;  
                            // collect edge
                            collected_edge.insert(std::pair<string, int>(this_contig_ID + "->" + tmp2_contig_ID + "_rev", 1));                                                               
                        }
                    }
                    */
                    double this_lod = calculate_lod(this_contig_GT, tmp2_contig_GT);
                    /*
                    if(this_lod>=40)
                    {
                        // select original: 	100052 -> 358538 [color=aquamarine, penwidth=1, arrowsize=1, label=45];
                        dotofp << "\t"   << this_contig_ID 
                               << " -> " << tmp2_contig_ID
                               << " "    << "[color=red, penwidth=1, arrowsize=1, label=" << this_lod << "];"    
                               << endl;
                        // collect edge
                        collected_edge.insert(std::pair<string, int>(this_contig_ID + "->" + tmp2_contig_ID, 1));                               
                    }
                    */
                    if(this_lod > contigBestLods.lod[0]) // best
                    {
                        contigBestLods.lod[0] = this_lod;
                        std::stringstream sstmp;
                        sstmp.str("");
                        sstmp << "\t"   << this_contig_ID 
                              << " -> " << tmp2_contig_ID
                              << " "    << "[color=red, penwidth=1, arrowsize=1, label=" << this_lod << "];";
                        contigBestLods.edge1  = sstmp.str();
                        // collect edge
                        collected_edge.insert(std::pair<string, int>(this_contig_ID + "->" + tmp2_contig_ID, 1));                                                       
                    }else
                    if(this_lod > contigBestLods.lod[1]) // second-best
                    {
                        //// contigBestLods.lod[1] = this_lod;
                        //// std::stringstream sstmp;
                        //// sstmp.str("");
                        //// sstmp << "\t"   << this_contig_ID 
                        ////      << " -> " << tmp2_contig_ID
                        ////      << " "    << "[color=red, penwidth=1, arrowsize=1, label=" << this_lod << "];";
                        //// contigBestLods.edge2  = sstmp.str(); 
                        //// // collect edge                        
                        //// collected_edge.insert(std::pair<string, int>(this_contig_ID + "->" + tmp2_contig_ID, 1));                                                                            
                    } else;
                }
            }        
            ctgitr2 ++;
        }
        //
        if(contigBestLods.lod[0] > 0) // best
        {
            dotofp << contigBestLods.edge1 << endl;
        }
        if(contigBestLods.lod[1] > 0) // second-best
        {
            //// dotofp << contigBestLods.edge2 << endl;         
        }
        //
        ctgitr ++;
    }
    // end of dot file
    dotofp << "}" << endl;
    dotofp.close();
    return true;
}
// check match score of genotypes of contigs 
bool get_match_score(string contig_GT1, string contig_GT2, int* score)
{
    if (contig_GT1.size() != contig_GT2.size()) 
    {
        cout << "   Error: sequences of contigs' genotypes are with different lengths. " << endl;
        return false;
    }
    *score = 0;
    for(int ii = 0; ii < contig_GT1.size(); ii ++)
    {
        if(contig_GT1[ii] == contig_GT2[ii])
        {
            *score += 1;
        }
    }
    return true;
}
// reverse GT of a contig
bool reverse_contig_GT(string contig_GT, string* contig_GT_rev)
{
    *contig_GT_rev = "";
    for(int ii = 0; ii < contig_GT.size(); ii ++)
    {
        if(contig_GT[ii] == 'P') 
        {
            *contig_GT_rev += "M"; 
        }else 
        if(contig_GT[ii] == 'M') 
        {
            *contig_GT_rev += "P"; 
        }else
        {
            *contig_GT_rev += contig_GT[ii]; 
        }
    }
    return true;
}
//
bool cluster_pollens(map<int, POLLEN>  pollens, 
                     map<int, POLLEN>* maCluster,
                     map<int, POLLEN>* paCluster,
                     map<int, POLLEN>* unCluster,
                     map<int, POLLEN>* pollens_updated)
{
    // idea: using initial haplotypes related to a contig as backbone, 
    //       1.color   pollen genotypes
    //       2.cluster pollen genotypes
    //       3.correct haplotypes
    assert(pollens.find(0) != pollens.end());
    string matGT = pollens[0].poSeq; // maternal backbone
    string patGT = pollens[0].poAll; // paternal backbone
    //
    // map<int, POLLEN> maCluster; // close to maternal GT
    // map<int, POLLEN> paCluster; // close to paternal GT
    // map<int, POLLEN> unCluster; // not informative 
    //
    map<int, POLLEN>::iterator pollenitr;
    map<int, POLLEN>::iterator pollenitr_end;
    pollenitr        = pollens.begin();
    pollenitr_end    = pollens.end();
    while(pollenitr != pollenitr_end)
    {
        int     pid  = (*pollenitr).first;
        POLLEN ptmp  = (*pollenitr).second;
        if(pid>0)
        {
            string        pollenGT = ((*pollenitr).second).poSeq;
            string        pattern  = "";                   
            unsigned long matscore = 0; 
            unsigned long patscore = 0;
            if( !get_similarty(matGT,
                               patGT,
                               pollenGT,
                               &pattern,
                               &matscore,
                               &patscore) )
            {
                return false;
            }
            // updating attribute pattern of pollen
            ((*pollenitr).second).poPat = pattern;
            // setting  attribute pattern of pollen
            ptmp.poPat = pattern;  
            //
            if(matscore  > patscore)
            {
                (*maCluster).insert(std::pair<int, POLLEN>(pid, ptmp));
            }
            else
            if(matscore  < patscore)
            {
                (*paCluster).insert(std::pair<int, POLLEN>(pid, ptmp));
            }else
            if(matscore == patscore)
            {
                (*unCluster).insert(std::pair<int, POLLEN>(pid, ptmp));
            }else ;
        }
        pollenitr ++;
    }
    // get pattern inserted in pollen data - new
    (*pollens_updated) = pollens;
    //
    return true;
}
bool get_similarty(string         matGT, 
                   string         patGT, 
                   string         pollenGT,
                   string*        pattern,                   
                   unsigned long* matscore, 
                   unsigned long* patscore)
{
    /*
      mat       : GATCCAGATGA
      pat       : AGCTTCATCAG
      pollenGT  : GANNCAGATGA
      =>        
      pattern =   MMuuMMMMMMM
      matscore= 9 (number of matched positions)
      patscore= 0 
    */
    if (matGT.size() != patGT.size()    || 
        matGT.size() != pollenGT.size() ||
        patGT.size() != pollenGT.size() ) 
    {
        cout << "   Error: sequences of pollen and haplotype initials are with different lengths. " << endl;
        return false;
    }
    //
    for(unsigned long ii=0; ii<matGT.size(); ii++)
    {
        if(matGT[ii]==pollenGT[ii])
        {
            (*matscore) += 1;
            (*pattern)  += "M";
        }else
        if(patGT[ii]==pollenGT[ii])
        {
            (*patscore) += 1;
            (*pattern)  += "P";
        }else
        {
            (*pattern)  += "u"; // undetermined/unclear
        }
    }
    //
    return true;
}
//
bool get_pollen_allele(string pollenfile, 
                       map<string, map<unsigned long, ALLELE> >* mkr)
{
    //
    ifstream ifp;
    ifp.open(pollenfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open pollen allele count file " << pollenfile << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;
        //contig	      pos     Consensus   totalcount  A	  C   G   T   -   N   ...........
        //tig00000013_pilon   14777   T           2           0   0   0   1   0   1   1.0     1.1
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 10) continue; // unexpected format
        string keyctg = lineinfo[0];
        // <ctg, <pos, {ref, alt}> >
        map<string, map<unsigned long, ALLELE> >::iterator ctgitr;
        ctgitr = (*mkr).find(keyctg);
        if(ctgitr == (*mkr).end() ) 
        {
            // pollen contig not in marker list
            continue;
        }
        else
        {
            //map<unsigned long, ALLELE> indMarker = (*ctgitr).second; 
            // directly update the pointers
            string                           keypos_str = lineinfo[1];
            unsigned long                        keypos = strtoul(keypos_str.c_str(), NULL, 0);
            map<unsigned long, ALLELE>::iterator positr = ( (*ctgitr).second ).find(keypos);
            if(positr == ( (*ctgitr).second ).end())
            {
                // pollen position not in marker list
                continue;
            }
            else
            {
                // find pollen genotype data for a marker
                string refb = (*positr).second.refb; // 
                string altb = (*positr).second.altb; // 
                int refpos = 4;
                if(!get_base_position(refb, &refpos)) return false;
                int altpos = 4;
                if(!get_base_position(altb, &altpos)) return false;
                unsigned long refcnt = strtoul(lineinfo[refpos].c_str(), NULL, 0);
                unsigned long altcnt = strtoul(lineinfo[altpos].c_str(), NULL, 0);
                // refb; altb; none of refb and altb found 'N'; both refb and altb found 'U'
                if(refcnt>0 && altcnt>0)
                {
                    (*positr).second.pollenb = "U";
                }
                else
                if(refcnt>0)
                {
                    (*positr).second.pollenb = refb;
                }
                else
                if(altcnt>0)
                {
                    (*positr).second.pollenb = altb;
                }
                else
                {
                    (*positr).second.pollenb = "N";
                }
            }
        }
    }
    ifp.close();
    return true;
}
//
bool get_base_position(string base, int* pos)
{
    // 4:A 5:C 6:G 7:T 8:- 9:N
    if(base.compare("A")==0) *pos = 4; else
    if(base.compare("C")==0) *pos = 5; else
    if(base.compare("G")==0) *pos = 6; else
    if(base.compare("T")==0) *pos = 7; else
    if(base.compare("-")==0) *pos = 8; else
    if(base.compare("N")==0) *pos = 9; else
    {
        cout << "   Error: unexpected base " << base << " in marker list" << endl;
        return false;
    }
    return true;
}
//
bool collect_contig_size(string sizefile, map<string, unsigned long>* contigsize)
{
    cout << "   Info: reading contig size info from " << sizefile << endl;
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
        cout << ctgid << "\t" << ctgsize << endl;
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
bool collect_marker_info(const char* file, 
                         map<string, map<unsigned long, ALLELE> >* mkr)
{
    // type of markers: snps(, 1-bp deletion and 1-bp insertion)
    cout << "   Info: reading marker info from " << file << endl;
    ifstream ifp;
    ifp.open(file);
    if(!ifp.good())
    {
        return false;
    }
    int num = 0;
    int raw = 0;
    int del = 0;
    int ins = 0;
    bool adjsnp = false; // SNPs directly next to each other
    int    nadj = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        /* e.g.,
	        ks	1	164393	C	G  -- snp
		ks	1	164413	AC	TA -- twp snps
		ks	1	460940  -	A  -- insertion
		ks	1	468923  T	-  -- deletion
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo[3].size() != lineinfo[4].size()) continue; // non-snps
        raw ++;
        
        if(lineinfo[3].size() > 1 && !adjsnp)
        {
            cout << "   Info: snp line will be separated, e.g, " << line << endl;
            cout << "       You will get a total number after reading all marker info. " << endl;
            adjsnp = true;
            nadj ++;
        }
        else
        if(lineinfo[3].size()>1)
        {
            nadj ++;
        }        
        //  split "double/triple/..." snps
        for(int snpi = 0; snpi < lineinfo[3].size(); snpi ++)
        {
            string ref = lineinfo[3].substr(snpi, 1);
            string alt = lineinfo[4].substr(snpi, 1);
            //
            unsigned long pos = strtoul(lineinfo[2].c_str(), NULL, 0);
            pos += snpi;
            std::stringstream sspos;
            sspos.str("");
            sspos << pos;
            //
            string        chrkey     = lineinfo[1];
            string        poskey_str = sspos.str();
            unsigned long poskey     = strtoul(poskey_str.c_str(), NULL, 0);
            string refb              = ref;
            string altb              = alt; // "GA", "G-" or "-A"
            string pollenb           = "N"; // unknown
            //
            map<string, map<unsigned long, ALLELE> >::iterator itr;
            map<string, map<unsigned long, ALLELE> >::iterator itr_end;
            itr     = (*mkr).find(chrkey);
            itr_end = (*mkr).end();
            
            if(itr == itr_end) 
            {
                ALLELE tmpallele;
                tmpallele.refb    = refb;
                tmpallele.altb    = altb;
                tmpallele.pollenb = pollenb;
                map<unsigned long, ALLELE> indMarker;
                indMarker.insert(std::pair<unsigned long, ALLELE>(poskey, tmpallele));
                // new contig
                (*mkr).insert(std::pair<string, map<unsigned long, ALLELE> >(chrkey, indMarker));
                num ++;
            } 
            else
            {
                ALLELE tmpallele;
                tmpallele.refb = refb;
                tmpallele.altb = altb;
                tmpallele.pollenb = pollenb;                
                // contig collected already with new marker position
                (*itr).second.insert(std::pair<unsigned long, ALLELE>(poskey, tmpallele));    
                num ++;                       
            }
            //
            if(alt=="-") del ++;
            if(ref=="-") ins ++;           
        }
    }
    ifp.close();
    cout << "   Info: recording " << num << " (of " << raw << ")" << " markers info done " << endl;  
    cout << "   Info: total number of \"neighboring\" snps: "     << nadj << endl;  
    cout << "   Info: total number of del markers: "              << del  << endl;
    cout << "   Info: total number of ins markers: "              << ins  << endl;
    if(num==0)
    {
        cout << "   Error: no markers found in marker file. " << endl;
        return false;
    }
    return true;
}
//
bool get_options(int                  argc, 
                 char*                argv[],
                 string*              filemark,
                 string*              filesize,
                 vector<string>*      filepoll,
                 bool*                correction,
                 double*              lodcutoff,
                 double*              minCOscore,
                 string*              outprefix)
{
    int  ic;
    bool markf   = false;
    bool pollf   = false;
    int  pollnum = 0; // number of pollen files given
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
            *filemark = (string)argv[ic];
            ifstream fp;
            fp.open((*filemark).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: marker file provided: "    << *filemark << endl;
            }
            else
            {
                cout << "   Error: cannot open marker file " << *filemark << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--corr") == 0)
        {
            *correction = true;
            cout << "   Info: correcting markers asked. " << endl;
        }
        else
        if(optstr.compare("--lod") == 0)
        {
            ic ++;
            *lodcutoff = atof(argv[ic]);
            cout << "   Info: lod cutoff provided to build conections: " << *lodcutoff << endl;
        } 
        else
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
        if(optstr.compare("--pollen") == 0)
        {
            ic ++;
            cout << "   Info: list of pollen variant files provided: "    << argv[ic] << endl;            
            ifstream fp;
            fp.open(argv[ic], ios::in);
            if(fp.good())
            {
                while(fp.good())
                {
                    string line("");
                    getline(fp, line);
                    if(line.size()==0 || line[0]=='#') continue;
                    pollnum ++;
                    ifstream ithfp;
                    ithfp.open(line.c_str(), ios::in);                    
                    if(ithfp.good())
                    {
                        (*filepoll).push_back(line);
                        ithfp.close();
                    }else
                    {
                        cout << "   Error: cannot open " << pollnum  << "th "  << "pollen file "  << line    << endl;
                    }                  
                }
                fp.close();
                cout << "   Info: No. pollen file collected/provided: " << (*filepoll).size() << "/" << pollnum << endl;
            }
            else
            {
                cout << "   Error: cannot open pollen-sample list file " << argv[ic] << endl;
                return false;
            }
        }
        else       
        if(optstr.compare("--size") == 0)
        {
            ic ++;
            ifstream fp;
            *filesize = (string)argv[ic];
            fp.open(argv[ic], ios::in);
            if(fp.good())
            {
                cout << "   Info: contig size file provided: "    << argv[ic] << endl;            
            }
            else
            {
                cout << "   Error: contig size file not found: "  << argv[ic] << endl;            
                return false;
            }
            fp.close();
        }        
        else
        if(optstr.compare("-o") == 0)
        {
            ic ++;
            *outprefix = (string)argv[ic];
            cout << "   Info: output files will be labeled with \""    << *outprefix << "\"" << endl;
        }        
        else
        {
            cout << "   Warning: option " << argv[ic] << " was not recognized and ignored."  << endl;
        }
        // next option
        ic ++;
    }
    // check necessary files
    return true;
}
