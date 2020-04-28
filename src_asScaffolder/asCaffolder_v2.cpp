/*  Given     
        a list of linkage map of ordered contigs -- from JoinMap              -- (X1:contigs of grouping+ordering),
        a list of .............. phased  contigs -- from JoinMap              -- (X2:contigs of grouping+phasing)        
        and         
        s2_genotype_contig_seq.txt               -- from asPollinator         -- (X3:all raw contigs     PM patterns) 
        s2_genotype_contig_seq_del_like.txt      -- from del_marker_genotyper -- (X4:all raw del-regions PM patterns)    
    Scaffold all contigs in linkage groups ( = complete linkage groups with PM-phasing and contig-ordering). 
    Note 1: contigs X1 <= a subset of X2 <= a subset of X3
    Note 2: contigs of X4 might be included in X1/2/3; 
            some=X5 of X4 might not =>X5 needs to be linkage-grouped, phased and ordered.    
    Note 3: contig key format: "contig\tleft" or "contig\tright" for X1/2/3; for X4, coordinates also considered.
    Note 4: on output, for each genetic map group, 
            P/M values assigned to del-like patterns can refer to different parents;
            as phasing can not be performed inter-groups. In the same genetic map group, P/M values should be 
            consistently naming the same parents.
    2020-02-11 13:28     updated: hap_low can be updated as hom_low according to cell support.
    Written by Hequan Sun, MPIPZ, Email: sunhequan@gmail.com
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
#include      <math.h>
#include "split_string.h"
using namespace std;
struct GNMAP
{
    // GeNetic MAP 
    string leri;    // left-end or right-end marker
    string parent;  // if in genetic map, who is the left marker
    string child;   // if in genetic map, who is the right marker
    double distance;// originally labeled genetic distance value 
    string group;   // linkage group
    string status;  // it is from JoinMap or updated with paired-end marker or PM pattern
};
// del-like marker
struct DELMARKER
{
    string        ctg;
    unsigned long sta;
    unsigned long end;       // end of del-like marker
    string        haphom;    // class: hap_high;hap_low;hom_low;hom_high
    string        pmpat;     // PM pattern for this del-like marker; if phase = 1, this becomes string of {M,U}
    string        gmap;      // which genetic map the marker would be assigned to
    string        status;    // indicate whether the pm pattern is changed
    long          score;     // score to put the marker into gmap
    string        insmarker; // before which snp-contig-marker, the del-like marker will be inserted.
};
double hap_low_to_hom_low_reset = 0.503; // 450/895=0.5027933
//
bool get_options(int                                          argc,
                 char*                                        argv[],
                 string*                                      fmap,
                 string*                                      fphase,
                 string*                                      fmarker,
                 string*                                      fmarker_del_like,
                 string*                                      outprefix, 
                 bool*                                        verbose);
bool get_group_files(string                                   listfiles, 
                 vector<string>*                              groupfiles,
                 string                                       gmapphasing);                   
bool get_contig_phase(string                                  jmfile, 
                 map<string, bool>*                           ctgJoinMapPhasing,
                 map<string, bool>*                           ctgJoinMapPhasing_all,
                 map<string, string>*                         ctgGroup,
                 int*                                         nflip); 
bool get_contig_genetic_map(string                            gmfile,
                 map<string, GNMAP>*                          gmContig,
                 vector<string>*                              gmContigOrder,
                 map<string, string>*                         unmappedContig);
bool get_contig_PMsequence(string                             mkrfile,
                 map<string, bool>                            ctgJoinMapPhasing_all,
                 map<string, string>*                         contigMarkerPM);
string retrieve_contigid(string                               jmcontigid); 
string flip_PMpat(string                                      PMpat);
bool match_phasegroup_mapgroup(map<string, vector<string> >   contig_genetic_map, 
                 map<string, map<string, bool> >              contig_phasing,
                 map<string, string>*                         matched_phase_map,
                 map<string, string>*                         matched_map_phase);
// complete 1: with snp-contig-markers
bool complete_genetic_map_snp_based(map<string, string>       contigMarkerPM, 
                 map<string, map<string, bool> >              contig_phasing,
                 map<string, string>                          ctgGroup,
                 map<string, vector<string> >                 contig_genetic_map, 
                 map<string, string>                          contig_gmap_group,                        
                 map<string, string>                          matched_phase_map,
                 map<string, string>                          matched_map_phase,
                 map<string, GNMAP>                           gmContig,
                 map<string, vector<string> >*                contig_genetic_map_updated,
                 map<string, bool>*                           ctgJoinMapPhasing_all,
                 map<string, string>*                         ctgGroup_updated,
                 map<string, string>*                         contig_gmap_group_updated,
                 map<string, string>*                         contigMarkerPM_updated);
// complete 2: with del-like-contig markers
bool complete_genetic_map_del_based(map<string, map<unsigned long, DELMARKER> > del_like_marker_PM, 
                 map<string, string>                          contigMarkerPM, 
                 map<string, map<string, bool> >              contig_phasing,
                 map<string, string>                          ctgGroup,
                 map<string, vector<string> >                 contig_genetic_map, 
                 map<string, string>                          contig_gmap_group,                        
                 map<string, string>                          matched_phase_map,
                 map<string, string>                          matched_map_phase,
                 map<string, GNMAP>                           gmContig,
                 map<string, vector<string> >*                contig_genetic_map_updated,
                 map<string, bool>*                           ctgJoinMapPhasing_all,
                 map<string, string>*                         ctgGroup_updated,
                 map<string, string>*                         contig_gmap_group_updated,
                 string                                       outprefix);              
long find_PM_pattern_match(string pm1,                        string pm2, 
                 bool                                         scale);     
string find_best_map_insertion_site(bool                      rely_on_paired,
                 string                                       this_raw_mkr, 
                 string                                       this_raw_pat,
                 vector<string>                               this_gm_contigs,
                 map<string, string>                          contigMarkerPM,
                 long*                                        highest_score);  // snp-contig-markers
string find_best_map_insertion_site_del_like(bool             rely_on_paired,
                 string                                       this_raw_mkr,
                 string                                       this_raw_pat, 
                 vector<string>                               this_gm_contigs,
                 map<string, string>                          contigMarkerPM,
                 long*                                        highest_score);  // del-like-contig-markers
double calculate_recomb_freq(string                           PMpat1, 
                 string                                       PMpat2);
// new after 20200102
bool get_del_like_PMsequence(string                           del_like_file,
                 map<string, map<unsigned long, DELMARKER> >* del_like_marker_PM);
bool output_del_like_PM_phasing(map<string, map<unsigned long, DELMARKER> > del_like_marker_PM)
{
    
};              
//
int main(int argc, char* argv[])
{
    if(argc < 11)
    {
        // g++ asCaffolder_v2.cpp split_string.cpp -O3 -o asCaffolder_v2
        cout << "\n   Function: scaffolding contigs using linkage-phased/ordered contig markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "asCaffolder_v2"
             << " --map linkage_map_list.txt"
             << " --phase ctg_phasing.txt"
             << " --marker s2_genotype_contig_seq.txt" 
             << " --marker-del-like s2_genotype_contig_seq_del_like.txt"
             << " --hap2hom 0.503"
             << " -o prefix_out"              << endl << endl;
        cout << "   Note: "                           << endl;
        cout << "      --map    is a file listing n files of X1: ordered contigs, n=linkgage group number. " << endl;        
        cout << "      --phase  is a file listing n files of X2: phased  contigs, n=linkgage group number. " << endl
             << "      --marker is a file listing X3: contig-PM patterns, e.g., s2_genotype_contig_seq.txt"  << endl
             << "               (X1 <= X2 <= X3.), given by asPollinator."                                   << endl
             << "      --marker-del-like is a file listing X4: contig-PM patterns" << endl
             << "               e.g., s2_genotype_contig_seq_del_like.txt by del_marker_genotyper. "         << endl
             << "      --hap2hom is a fraction of cells support the region as del or not - for 895 cells, "  << endl
             << "                we require less than half cells; if more than that, a hap_low would be reset"<< endl
             << "                as hom_low which influence pacbio read separation to 1 or 2 haplotypes. "   << endl               
             << endl;
        return 1;
    }
    double startT= clock();
    // step 0. initialize from the cmd line
    cout << "Initialize input from cmd line: " << endl << endl;
    string fmap             = "";
    string fphase           = "";    
    string fmarker          = "";
    string fmarker_del_like = "";
    string outprefix        = "fun_";
    bool   verbose          = false;
    if(!get_options(argc,
                    argv,
                    &fmap,
                    &fphase,                    
                    &fmarker,
                    &fmarker_del_like,
                    &outprefix, 
                    &verbose) )
    {
        return 1;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// step 1. get contig map with limited contigs - backbone of scaffolding
    //////////////////////////////////         with 1st completing of the genetic map: if the paired marker is mapped 
    //////////////////////////////////              insert the unmapped marker next to its paired marker /////////////// 
    cout << " Step 1: reading genetic map info on contig markers..." << endl;
    map<string, vector<string> > contig_genetic_map;// map<this_map_file, vector<contig-id-in-map> > // file-grouped map
    map<string, GNMAP> gmContig;    
    vector<string> groupmap;
    if(!get_group_files(fmap, &groupmap, "genetic map"))
    {
        return 1;
    }
    vector<string>::iterator mapfitr;
    vector<string>::iterator mapfitr_end;
    mapfitr     = groupmap.begin();
    mapfitr_end = groupmap.end();
    while(mapfitr != mapfitr_end)
    {
        string this_gmfile           = *mapfitr;  
        vector<string> fileinfo      = split_string(this_gmfile, '/');
        string this_genetic_map_file = fileinfo[fileinfo.size()-1];  // only file name without path info              
        cout << endl;
        cout << "   Info: reading genetic map info from " << this_gmfile << "..."    << endl;        
        vector<string>       gmContigOrder;
        map<string, string>  unmappedContig;                   
        if(!get_contig_genetic_map(this_gmfile,
                                   &gmContig,
                                   &gmContigOrder,
                                   &unmappedContig))
        {
            cout << "   Error: failed in reading " << this_gmfile << endl;
            return 1;
        }
        cout << "   Info: reading genetic map info from " << this_gmfile << " done. " << endl;       
        //
        contig_genetic_map.insert(std::pair<string, vector<string> >(this_genetic_map_file, gmContigOrder));      
        //
        mapfitr ++;
    }
    cout << " Step 1: reading genetic map info on contig markers done."    << endl  << endl;       
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// step 2. get contig phasing
    cout << " Step 2: reading phasing info on contig markers in linkage groups... " << endl;
    vector<string> groupphase;      
    if(!get_group_files(fphase, &groupphase, "phasing/grouping"))
    {
        return 1;
    }
    map<string, map<string, bool> > contig_phasing;        // map<this_linkage_file, map<contig_id, bool_flip > >
    map<string, bool>               ctgJoinMapPhasing_all; // map<contig_id, bool-flip>
    map<string, string>             ctgGroup;              // map<contig_id, this_linksage_file>
    int                             nflip_all = 0;         // total number of PM patterns of contig markers need flipping
    int                             nflip     = 0;         // group-number of PM patterns of contig markers need flipping    
    //
    cout << endl;
    vector<string>::iterator gpitr;
    vector<string>::iterator gpitr_end;
    gpitr     = groupphase.begin();
    gpitr_end = groupphase.end();
    while(gpitr != gpitr_end)
    {
        string this_jmfile = *gpitr;
        vector<string> fileinfo  = split_string(this_jmfile, '/');
        string this_linkage_file = fileinfo[fileinfo.size()-1];        
        cout << "   Info: read linkage-phasing info from " << fileinfo[fileinfo.size()-1] << endl;    
        map<string, bool> ctgJoinMapPhasing;
        if(!get_contig_phase(this_jmfile, 
                             &ctgJoinMapPhasing,
                             &ctgJoinMapPhasing_all,
                             &ctgGroup,
                             &nflip))
        {
            cout << "   Error: failed in reading " << this_jmfile << endl;
            return 1;
        }
        nflip_all += nflip;
        contig_phasing.insert(std::pair<string, map<string, bool> >(this_linkage_file, ctgJoinMapPhasing));
        cout << "         phasing info on " << ctgJoinMapPhasing.size() << " contig markers collected. " << endl;
        //
        gpitr ++;
    }
    cout << "   Info: overall, phasing info on " << ctgJoinMapPhasing_all.size() << " contig markers collected." << endl;
    cout << "   Info: overall, PM pattern of "   << nflip_all << " contig markers need flipping. "         << endl;
    cout << " Step 2: reading phasing info on contig markers in linkage groups done. "               << endl << endl;   
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    /////////////////////////////////////// step 3: match phased group and genetic map group
    cout << " Step 3: matching the phased group of contigs with groups of contigs in genetic map... "       << endl;    
    map<string, string> matched_phase_map; // <phasefile,  gmapfile> 
    map<string, string> matched_map_phase; // <gmapfile,  phasefile>
    if(!match_phasegroup_mapgroup(contig_genetic_map, 
                                  contig_phasing,
                                  &matched_phase_map,
                                  &matched_map_phase))
    {
        cout << "   Error: cannot match phased contigs with those in genetic map. " << endl;
        return 1;
    }    
    if(1)
    {
        cout << "   Info: phased group should be more complete than mapped regarding number of markers. "     << endl;    
        map<string, string>::iterator checkitr;
        map<string, string>::iterator checkitr_end;
        checkitr     = matched_phase_map.begin();
        checkitr_end = matched_phase_map.end();
        while(checkitr != checkitr_end)
        {
            cout << "   Info: phased-group: "                 << (*checkitr).first  
                 << " <=paired with=> mapped/ordered-group: " << (*checkitr).second << endl;
            checkitr ++;
        }
    }
    cout << " Step 3: matching the phased group of contigs with groups of contigs in genetic map done. "<<endl<< endl;    
         
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// step 4. get contig PM patterns from asPollinator
    cout << " Step 4: reading PM patterns on contig markers from asPollinator... "            << endl;    
    // note: contigMarkerPM would collect PM patterns flipped based on phase info from JoinMap
    map<string, string> contigMarkerPM; //  <contig\tleft/right, PMpattern>
    if(!get_contig_PMsequence(fmarker, ctgJoinMapPhasing_all, &contigMarkerPM))
    {
        cout << "   Error: failed in reading file " << fmarker << endl;
        return 1;
    }
    cout << " Step 4: reading PM patterns on contig markers from asPollinator done. " << endl << endl;    
    // step 4.5 finding intermediate info about contig id and group of genetic map
    map<string, string> contig_gmap_group;   // <contig\tleft/right, this_genetic_map_file>
    map<string, vector<string> >::iterator mfitr;
    map<string, vector<string> >::iterator mfitr_end;
    mfitr     = contig_genetic_map.begin();
    mfitr_end = contig_genetic_map.end();
    while(mfitr != mfitr_end)
    {
        string this_gm_file = (*mfitr).first;
        vector<string> this_gm_contigs = (*mfitr).second;
        vector<string>::iterator gmcitr;
        vector<string>::iterator gmcitr_end;
        gmcitr     = this_gm_contigs.begin();
        gmcitr_end = this_gm_contigs.end();
        while(gmcitr != gmcitr_end)
        {
            string this_contig = *gmcitr;
            contig_gmap_group.insert(std::pair<string, string>(this_contig, this_gm_file));
            gmcitr ++;
        }
        mfitr ++;
    }
    cout << "   Info: " << contig_gmap_group.size() << " contig markers in genetic map - reminder. " << endl   << endl;
    cout << " Step 4.x: reading PM patterns on del-like contig-regional markers from del_marker_genotyper... " << endl;
    map<string, map<unsigned long, DELMARKER> > del_like_marker_PM;
    if(!get_del_like_PMsequence(fmarker_del_like, &del_like_marker_PM))
    {
        cout << "   Error: failed in reading file " << fmarker_del_like << endl;
        return 1;
    }
    cout << "   Info: such del-like markers need to be inserted into the genetic map based on PM patterns later: \n"
         << "         if the related contig is in genetic map, "                   << endl
         << "              if it is hap-like region, only need phasing; "          << endl
         << "              if it is hom-like region, needing no action (=>pacbio reads go to both genotypes); " << endl
         << "         else,  "                                                     << endl
         << "              if it is hap-like reigon, need grouping and phasing;"   << endl
         << "              if it is hom-like reigon, only need grouping (=>pacbio reads go to both genotypes)." << endl;
    cout << " Step 4.x: reading PM patterns on del-like contig-regional markers from del_marker_genotyper done. "
         << endl
         << endl;
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /////////////////////////////////////// step 5. 2nd completing on genetic map    
    // 2nd completing of the genetic map according PM pattern similarity of snp-defined contig markers.
    cout << " Step 5: completing genetic map with PM patterns of contig markers from asPollinator... "     << endl; 
    cout << "         Note: for a missing marker x:X(PM-pattern), and candidate insertion sites y:Y and z:Z: "
         << endl
         << "               ins score  = 0;"                                                               << endl           
         << "               ins score += num-matched-pos(X,Y) / num-both-XYpos-are-P-or-M (X,Y) * X-size;" << endl       
         << "               ins score += num-matched-pos(X,Z) / num-both-XZpos-are-P-or-M (X,Z) * X-size." << endl
         << endl;    
    map<string, vector<string> >   contig_genetic_map_updated;    
    map<string, string>            ctgGroup_updated;
    map<string, string>            contig_gmap_group_updated;
    map<string, string>            contigMarkerPM_updated;
    if(!complete_genetic_map_snp_based(contigMarkerPM, 
                                       contig_phasing,
                                       ctgGroup,
                                       contig_genetic_map, 
                                       contig_gmap_group,
                                       matched_phase_map,
                                       matched_map_phase,                                       
                                       gmContig,                             
                                       &contig_genetic_map_updated,
                                       &ctgJoinMapPhasing_all,
                                       &ctgGroup_updated,
                                       &contig_gmap_group_updated,
                                       &contigMarkerPM_updated))
    {
        cout << "   Error: failed in completing the genetic map. " << endl;
        return 1;
    }    
    cout << " Step 5: completing genetic map with PM patterns on contig markers from asPollinator done. " 
         << endl << endl;    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /////////////////////////////////////// step 6. 3rd completing on genetic map with del-like markers   
    // 3rd completing of the genetic map according PM pattern similarity of del-like contig-regional markers.
    cout << " Step 6: completing genetic map with PM patterns of del-like contig markers from del_marker_genotyper... "
         << endl;     
    ctgGroup.clear();
    contig_genetic_map.clear();
    contig_gmap_group.clear();
    contigMarkerPM.clear();
    ctgGroup           = ctgGroup_updated;
    contig_genetic_map = contig_genetic_map_updated;
    contig_gmap_group  = contig_gmap_group_updated;
    contigMarkerPM     = contigMarkerPM_updated;    
    contig_genetic_map_updated.clear();    
    contig_gmap_group_updated.clear();    
    if(!complete_genetic_map_del_based(del_like_marker_PM, 
                                       contigMarkerPM,     
                                       contig_phasing,
                                       ctgGroup,
                                       contig_genetic_map, 
                                       contig_gmap_group,
                                       matched_phase_map,
                                       matched_map_phase,                                       
                                       gmContig,                             
                                       &contig_genetic_map_updated,
                                       &ctgJoinMapPhasing_all,
                                       &ctgGroup_updated,
                                       &contig_gmap_group_updated,
                                       outprefix                                       
                                      ))
    {
        cout << "   Error: failed in completing the genetic map. " << endl;
        return 1;
    } 
    cout << " Step 6: completing genetic map with PM patterns of del-like contig markers from del_marker_genotyper done."
         << endl << endl;        
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
// complete 3: complete genetic map with del-like contigs 
//             (if the contig is not in the genetic map corrected with snp-contig-marker)
bool complete_genetic_map_del_based(map<string, map<unsigned long, DELMARKER> >  del_like_marker_PM, 
                                    map<string, string>                          contigMarkerPM, 
                                    map<string, map<string, bool> >              contig_phasing,
                                    map<string, string>                          ctgGroup,
                                    map<string, vector<string> >                 contig_genetic_map, 
                                    map<string, string>                          contig_gmap_group,                        
                                    map<string, string>                          matched_phase_map,
                                    map<string, string>                          matched_map_phase,
                                    map<string, GNMAP>                           gmContig,
                                    map<string, vector<string> >*                contig_genetic_map_updated,
                                    map<string, bool>*                           ctgJoinMapPhasing_all,
                                    map<string, string>*                         ctgGroup_updated,
                                    map<string, string>*                         contig_gmap_group_updated,
                                    string                                       outprefix)
{
    /* 
    allraw del-like contig regional markers    - del_like_marker_PM   : <contig, <start, {ctg, sta, end, haphom, pmpat, group, status}> >
    -- 
    -- grouping related (from JoinMap)
    phased contig-left/right markers in groups - contig_phasing       : <this_linkage_file , <contig\tleft/right, bool_flip> > 
    phased contig-left/right markers           - ctgGroup             : <contig\tleft/right, this_linksage_file>
    UPDATED contig phasing information         - ctgJoinMapPhasing_all: <contig\tleft/right, bool_flip>    
    -- 
    -- ordering related (from JoinMap)
    mapped contig-left/right markers           - contig_genetic_map   : <this_gmap_file    , <contig\tleft/right> >
    mapped contig-left/right markers in maps   - contig_gmap_group    : <contig\tleft/right, this_gmap_file> 
    -- 
    -- 
    matching info of phased file and map file  - matched_phase_map    : <this_linkage_file , this_gmap_file>  
    -- 
    -- updated genetic map with snp-contig-markers
    UPDATED mapped contig-left/right markers   - contig_genetic_map_updated: <this_gmap_file    , <contig\tleft/right> >  
    */ 
    // initialize the genetic map with mapped contig markers
    // ctgJoinMapPhasing_all: this includes initial phased contigs; to be updated by integrating ungrouped_ctgPhasing    
    ////(*contig_genetic_map_updated).insert(contig_genetic_map.begin(), contig_genetic_map.end());
    ////(*ctgGroup_updated).insert(ctgGroup.begin(), ctgGroup.end());
    ////(*contig_gmap_group_updated).insert(contig_gmap_group.begin(), contig_gmap_group.end());
    map<string, bool> ungrouped_ctgPhasing;  
    map<string, string> contigMarkerPM_updated;
    ////contigMarkerPM_updated.insert(contigMarkerPM.begin(), contigMarkerPM.end());
    //
    map<string, string> to_be_inserted_del_contig_id;    // <del-contig, best_gmapfile+"#"+best_snp-contig>
    map<string, string> to_be_inserted_del_contig_pmpat; // <del-contig, best_representative_pmpat>        
    // prepare output file for collecting phased info 
    string ophasedfile = outprefix + "_s2_genotype_contig_seq_del_like.txt";
    ofstream ofp;
    ofp.open(ophasedfile.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file to collected phased del-like markers: " << ophasedfile << endl;
        return false;
    }
    //
    map<string, map<unsigned long, DELMARKER> >::iterator delctgitr;
    map<string, map<unsigned long, DELMARKER> >::iterator delctgitr_end;
    delctgitr     = del_like_marker_PM.begin();
    delctgitr_end = del_like_marker_PM.end();
    unsigned long case1_hom_cnt = 0;
    unsigned long case1_hap_cnt = 0;
    unsigned long case1_bad_cnt = 0; // no P value in pmpat
    unsigned long case2_hom_cnt = 0;
    unsigned long case2_hap_cnt = 0;
    unsigned long case3_hom_cnt = 0;
    unsigned long case3_hap_cnt = 0;    
    unsigned long caseall_total = 0; 
    unsigned long not_in_gmap_not_in_phase = 0;
    unsigned long not_in_gmap_but_in_phase = 0;
    unsigned long in_gmap_thus_in_phase    = 0;
    map<int, map<string, int> > ctg_case123;
    map<string, int> tmpcnt;
    ctg_case123.insert(std::pair<int, map<string, int> >(1, tmpcnt)); // case 1 contigs
    ctg_case123.insert(std::pair<int, map<string, int> >(2, tmpcnt)); // case 2 contigs
    ctg_case123.insert(std::pair<int, map<string, int> >(3, tmpcnt)); // case 3 contigs
    //
    while(delctgitr != delctgitr_end)
    {
        string this_ctg = (*delctgitr).first;
        string this_left_id         = this_ctg + "\tleft"; // building ids like ctg\tleft/right in asPollinator format
        string this_right_id        = this_ctg + "\tright";           
        string this_linkage_file = "";                
        if(ctgGroup.find(this_left_id)  != ctgGroup.end())
        {
            this_linkage_file = ctgGroup[this_left_id];
        }else
        if(ctgGroup.find(this_right_id)  != ctgGroup.end())
        {
            this_linkage_file = ctgGroup[this_right_id];
        }else ;
        // a list of del-like for current contig: 
        //    all del-like from this contig should go to the same genetic map/ordering. 
        //    if all hom-like, we cannot do anything; if at least one hap-like, it can be groupped.           
        map<unsigned long, DELMARKER> this_del_like = (*delctgitr).second;     
        map<unsigned long, DELMARKER>::iterator ditr;
        map<unsigned long, DELMARKER>::iterator ditr_end;
        ditr     = this_del_like.begin();
        ditr_end = this_del_like.end();
        map<string, int>  gmap_count; // how many cases go to which genetic map.
        map<string, int>  gmap_score; // maximum score of the marker going to a genetic map.
        map<string, map<string, int> >     mark_score;       // <gmapfile, <contig_known_marker, insertion_score> >
        map<string, map<string, string> >  mark_score_pmpat; // <gmapfile, <contig_known_marker+"\t"+del_marker, insertion_del_marker_pmpat> >
        ////map<string, unsigned long> gmap_start; // sta with maximum score of the marker going to a genetic map.
        while(ditr != ditr_end)
        {        
            //
            caseall_total ++;            
            //
            DELMARKER this_node         = (*ditr).second;      // {ctg, sta, end, haphom, pmpat, gmap, status, insmarker}    
            unsigned long this_del_sta  = this_node.sta;
            //
            string this_raw_pat         = this_node.pmpat;
            string this_raw_mkr         = this_left_id;        // left/right have no effect as rely_on_paired=false later
            string this_raw_pat_flipped = flip_PMpat(this_raw_pat);          
            // special contigs
            if(this_ctg.compare("tig00003688_pilon")==0 || this_ctg.compare("tig00003826_pilon")==0)
            {
                this_left_id  = this_ctg + "sA\tleft";
                this_right_id = this_ctg + "sA\tright";
            }
            //
            std::stringstream ss;                             // building id in del_marker_genotyper format
            ss.str("");
            ss << this_ctg << "\t" << this_node.sta << "\t" << this_node.end;
            // level 1: check if current marker in genetic map or not = ordering checking.
            if(contig_gmap_group.find(this_left_id)  == contig_gmap_group.end() && 
               contig_gmap_group.find(this_right_id) == contig_gmap_group.end())
            {
                cout << "       : case 1or2 - contig of this del-like marker \"" << ss.str() << "\" NOT in genetic map: " << endl;
                // level 2: check if current marker in phased group or not = phasing checking
                if(ctgGroup.find(this_left_id)  == ctgGroup.end() &&
                   ctgGroup.find(this_right_id) == ctgGroup.end())
                {
                    not_in_gmap_not_in_phase ++;    
                    if(ctg_case123[1].find(this_ctg) == ctg_case123[1].end())
                    {
                       ctg_case123[1].insert(std::pair<string, int>(this_ctg, 1));
                    }              
                    // case 1: not in ordering + not in grouping/phasing
                    cout << "       : case 1 - also NOT in phased groups, "                     << endl;   
                    if(this_node.haphom.find("hom") != std::string::npos)
                    {
                        // case 1.1: hom 
                        cout << "                  ";
                        cout << "but its depth>=hom => no further genetic map completing can be done with it; "
                             << "pacbio reads overlapping such a region would go to unmapped group. " << endl << endl;
                        this_del_like[this_del_sta].gmap   = "gmap.file.unavailable";
                        this_del_like[this_del_sta].status = "unphased.hom.case1.1";
                        this_del_like[this_del_sta].score  = 0;
                        if(gmap_count.find((string)"gmap.file.unavailable") == gmap_count.end())
                        {
                            gmap_count.insert(std::pair<string, int>("gmap.file.unavailable", 1));
                            gmap_score.insert(std::pair<string, int>("gmap.file.unavailable", 0));
                            ////gmap_start.insert(std::pair<string, unsigned long>("gmap.file.unavailable", this_del_sta));
                        }else
                        {
                            gmap_count["gmap.file.unavailable"]    += 1;
                            // no need to update gmap_score in this case as it is always 0.
                            ////size_t n = std::count(this_raw_mkr.begin(), this_raw_mkr.end(), 'P');
                            ////unsigned long last_sta = gmap_start["gmap.file.unavailable"];
                            ////size_t m = std::count(this_del_like[last_sta].pmpat.begin(), this_del_like[last_sta].pmpat.end(), 'P');
                            ////if(n > m)
                            ////{
                            ////    gmap_start.insert(std::pair<string, unsigned long>("gmap.file.unavailable", this_del_sta));                                
                            ////}
                        }
                        case1_hom_cnt ++;
                    }else
                    {
                        // case 1.2: hap
                        cout << "                  ";                        
                        cout << "but it is hap => phasing it according PM similarity "
                             << "to contigs in all genetic maps to decide genotype of pacbio reads. " << endl;
                        // case 1.2.1
                        if(this_raw_pat.find("P") == std::string::npos)
                        {
                            cout << "                  "; 
                            cout << "there is NO read info for all pollen nuclei; "
                                 << "pacbio reads overlapping such a region would go to unmapped group. " <<endl<< endl; 
                            this_del_like[this_del_sta].gmap   = "gmap.file.unavailable";
                            this_del_like[this_del_sta].status = "unphased.hap.case1.2.1";  
                            this_del_like[this_del_sta].score  = 0;                            
                            if(gmap_count.find((string)"gmap.file.unavailable") == gmap_count.end())
                            {
                                gmap_count.insert(std::pair<string, int>("gmap.file.unavailable", 1));
                                gmap_score.insert(std::pair<string, int>("gmap.file.unavailable", 0));  
                                ////gmap_start.insert(std::pair<string, unsigned long>("gmap.file.unavailable", this_del_sta));                                                        
                            }else
                            {
                                gmap_count["gmap.file.unavailable"]    += 1;
                                // no need to update gmap_score in this case as it is always 0.    
                                ////size_t n = std::count(this_raw_mkr.begin(), this_raw_mkr.end(), 'P');
                                ////unsigned long last_sta = gmap_start["gmap.file.unavailable"];
                                ////size_t m = std::count(this_del_like[last_sta].pmpat.begin(), this_del_like[last_sta].pmpat.end(), 'P');
                                ////if(n > m)
                                ////{
                                ////    gmap_start.insert(std::pair<string, unsigned long>("gmap.file.unavailable", this_del_sta));                                
                                ////}                                                            
                            }
                            case1_bad_cnt ++;                                                            
                        }else
                        {                             
                            // find insertion site of current marker in the genetic map
                            // target contigs defined with snp-markers in all genetic maps: contig_genetic_map
                            map<string, vector<string> >::iterator gmfitr;
                            map<string, vector<string> >::iterator gmfitr_end;
                            gmfitr     = contig_genetic_map.begin();
                            gmfitr_end = contig_genetic_map.end();
                            long highest_score_overall = 0;
                            string insertion_site_ctg_overall   = "";
                            string this_gmap_file_overall       = "";
                            bool phase_assign                   = false;                            
                            while(gmfitr != gmfitr_end)
                            {
                                cout << "       : investigating genetic map " << (*gmfitr).first << endl;
                                string         this_gmap_file   = (*gmfitr).first;
                                vector<string> this_gm_contigs  = (*gmfitr).second;
                                // raw
                                long highest_score;
                                bool rely_on_paired = false;
                                string insertion_site_ctg = find_best_map_insertion_site_del_like(rely_on_paired,
                                                                                         this_raw_mkr,
                                                                                         this_raw_pat, 
                                                                                         this_gm_contigs,
                                                                                         contigMarkerPM,
                                                                                         &highest_score);
                                cout << "       : score with P pattern = " << highest_score << endl;                                                                                         
                                if(highest_score == highest_score_overall)
                                {
                                    cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                                         << highest_score_overall << endl;
                                }                                                                             
                                if(highest_score > highest_score_overall)
                                {
                                    highest_score_overall      = highest_score;
                                    insertion_site_ctg_overall = insertion_site_ctg;
                                    this_gmap_file_overall     = this_gmap_file;
                                    phase_assign               = false;
                                }
                                // flipped                                                          
                                long highest_score_flipped;
                                rely_on_paired = false;                    
                                string insertion_site_ctg_flipped = find_best_map_insertion_site_del_like(rely_on_paired,
                                                                                         this_raw_mkr,
                                                                                         this_raw_pat_flipped, 
                                                                                         this_gm_contigs,
                                                                                         contigMarkerPM,
                                                                                         &highest_score_flipped);
                                cout << "       : score with M pattern = " << highest_score_flipped << " (flipped) "<< endl;
                                if(highest_score_flipped == highest_score_overall)
                                {
                                    cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                                         << highest_score_overall << endl;
                                }                                                                             
                                if(highest_score_flipped > highest_score_overall)
                                {
                                    highest_score_overall      = highest_score_flipped;
                                    insertion_site_ctg_overall = insertion_site_ctg_flipped;
                                    this_gmap_file_overall     = this_gmap_file;          
                                    phase_assign               = true;              
                                }                                                          
                                //
                                gmfitr ++;
                            }      
                            //
                            vector<string>::iterator ins_itr = std::find(contig_genetic_map[this_gmap_file_overall].begin(),
                                                                         contig_genetic_map[this_gmap_file_overall].end(),
                                                                         insertion_site_ctg_overall);              
                            // note del-like can be multiple for the same contig, so "sta\tend" as the label. 
                            if(phase_assign == true)
                            {                         
                                // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat_flipped)); 
                                cout << this_raw_pat_flipped     << "\t" << ss.str() << endl
                                     << contigMarkerPM[*ins_itr] << "\t" << *ins_itr << endl;
                                this_raw_pat = this_raw_pat_flipped;
                            }else
                            {
                                // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat));
                                cout << this_raw_pat             << "\t" << ss.str() << endl
                                     << contigMarkerPM[*ins_itr] << "\t" << *ins_itr << endl;                                
                            }
                            //
                            // update genetic map and related information
                            cout << "       : inserting before " << *ins_itr
                                 << " with score "               << highest_score_overall 
                                 << " to  genetic map group "    << this_gmap_file_overall 
                                 << " phase={"                   << phase_assign << "}"
                                 << endl;                             
                            // contig_genetic_map[this_gmap_file_overall].insert(ins_itr, ss.str());                             
                            // new phase status for ungrouped contigs
                            // ungrouped_ctgPhasing.insert(std::pair<string, bool>(this_raw_mkr, phase_assign));    
                            // update
                            // contig_gmap_group.insert(std::pair<string, string>(ss.str(), 
                               //                                                this_gmap_file_overall));
                            // update
                            // assert(matched_map_phase.find(this_gmap_file_overall) != matched_map_phase.end());
                            cout << "       : related PM phasing file: " << matched_map_phase[this_gmap_file_overall] 
                                 << ". " << endl << endl;                 
                            // (*ctgGroup_updated).insert(std::pair<string, string>(this_raw_mkr, 
                            //                                                     matched_map_phase[this_gmap_file_overall] ));
                            // update inserted linkage group   
                            this_del_like[this_del_sta].gmap      = this_gmap_file_overall;
                            this_del_like[this_del_sta].status    = "phased.hap.case1.2.2";
                            this_del_like[this_del_sta].pmpat     = this_raw_pat; // might have become "M" pattern
                            this_del_like[this_del_sta].insmarker = *ins_itr;
                            this_del_like[this_del_sta].score     = highest_score_overall;   
                            //     
                            if(mark_score.find(this_gmap_file_overall) == mark_score.end())
                            {
                                map<string, int> tmpscore;
                                tmpscore.insert(std::pair<string, int>(*ins_itr, highest_score_overall));
                                mark_score.insert(std::pair<string, map<string, int> >(this_gmap_file_overall, 
                                                                                       tmpscore));
                                map<string, string> tmppmpat;
                                tmppmpat.insert(std::pair<string, string>(*ins_itr + "\t" + this_ctg, this_raw_pat));
                                mark_score_pmpat.insert(std::pair<string, map<string, string> >(this_gmap_file_overall, 
                                                                                                tmppmpat));                                                      
                            }else
                            {
                                if(mark_score[this_gmap_file_overall].find(*ins_itr) == mark_score[this_gmap_file_overall].end())
                                {
                                    mark_score[this_gmap_file_overall].insert(std::pair<string, int>(*ins_itr, highest_score_overall));
                                    mark_score_pmpat[this_gmap_file_overall].insert(std::pair<string, string>(*ins_itr + "\t" + this_ctg, this_raw_pat));
                                }else
                                {
                                    if(highest_score_overall > mark_score[this_gmap_file_overall][*ins_itr])
                                    {
                                        mark_score[this_gmap_file_overall][*ins_itr] = highest_score_overall;
                                        mark_score_pmpat[this_gmap_file_overall][*ins_itr + "\t" + this_ctg] = this_raw_pat;
                                    }
                                }
                            }
                            //
                            if(gmap_count.find(this_gmap_file_overall) == gmap_count.end())
                            {
                                gmap_count.insert(std::pair<string, int>(this_gmap_file_overall, 1));
                                gmap_score.insert(std::pair<string, int>(this_gmap_file_overall, highest_score_overall));  
                                ////gmap_start.insert(std::pair<string, unsigned long>(this_gmap_file_overall, this_del_sta));
                            }else
                            {
                                gmap_count[this_gmap_file_overall] += 1;
                                if(gmap_score[this_gmap_file_overall] < highest_score_overall)
                                {
                                    gmap_score[this_gmap_file_overall] = highest_score_overall;                                
                                    ////gmap_start[this_gmap_file_overall] = this_del_sta;
                                }
                            }
                            //
                            case1_hap_cnt ++;
                        }
                    }
                }else
                {
                    not_in_gmap_but_in_phase ++;
                    if(ctg_case123[2].find(this_ctg) == ctg_case123[2].end())
                    {
                       ctg_case123[2].insert(std::pair<string, int>(this_ctg, 1));
                    }                     
                    // case 2: not in genetic map/ordering but in grouping/phasing
                    string this_linkage_file = "";
                    if(ctgGroup.find(this_left_id)  != ctgGroup.end())
                    {
                        this_linkage_file = ctgGroup[this_left_id];
                    }else
                    if(ctgGroup.find(this_right_id)  != ctgGroup.end())
                    {
                        this_linkage_file = ctgGroup[this_right_id];
                    }else ;
                    
                    string this_gmap_file          = matched_phase_map[this_linkage_file];
                    vector<string> this_gm_contigs = contig_genetic_map[this_gmap_file];                
                    cout << "       : case 2- but IN phased group "    << this_linkage_file 
                         << "; matching genetic map file "             << this_gmap_file
                         << " with " << this_gm_contigs.size()         << " ctg-markers. " 
                         << endl;  
                    if(this_node.haphom.find("hom") != std::string::npos)
                    {
                        // case 2.1: hom 
                        cout << "                  ";                                                
                        cout << "but its depth>=hom => no further genetic map completing can be done with it; "
                             << "pacbio reads overlapping such a region would go to both genotypes of group "
                             << this_linkage_file  
                             << endl;
                        this_del_like[this_del_sta].gmap   = this_gmap_file;
                        this_del_like[this_del_sta].status = "unphased.hom.case2.1";  
                        this_del_like[this_del_sta].score  = 2*this_raw_pat.size();
                        //
                        if(gmap_count.find(this_gmap_file) == gmap_count.end())
                        {
                            gmap_count.insert(std::pair<string, int>(this_gmap_file, 1));
                            gmap_score.insert(std::pair<string, int>(this_gmap_file, 2*this_raw_pat.size()));   
                            ////gmap_start.insert(std::pair<string, unsigned long>(this_gmap_file, this_del_sta));                                                     
                        }else
                        {
                            gmap_count[this_gmap_file] += 1;
                            if(gmap_score[this_gmap_file] < 2*this_raw_pat.size())
                            {
                                gmap_score[this_gmap_file] = 2*this_raw_pat.size();
                                ////gmap_start[this_gmap_file] = this_del_sta;                                
                            }
                        }
                        case2_hom_cnt ++;                                                   
                    }else
                    {
                        // case 2.2: hap
                        cout << "                  ";                        
                        cout << "and it is hap => phasing it according PM similarity "
                             << "to " << this_ctg << " in group "             
                             << this_linkage_file
                             << " to decide genotype of pacbio reads. " 
                             << endl;
                        cout << "                  note: if " << this_ctg << " does not have enough effective pollen; "
                             << "using other PM of other contigs in the same group. "
                             << endl;
                        //
                        cout << "       : investigating genetic map " << this_gmap_file << endl;
                        // find insertion site of current marker in the genetic map
                        long highest_score;
                        bool rely_on_paired = false;
                        long highest_score_overall = 0;
                        string insertion_site_ctg_overall   = "";
                        bool phase_assign                   = false;                            
                        string insertion_site_ctg = find_best_map_insertion_site_del_like(rely_on_paired,
                                                                                         this_raw_mkr,
                                                                                         this_raw_pat, 
                                                                                         this_gm_contigs,
                                                                                         contigMarkerPM,
                                                                                         &highest_score);
                        cout << "       : score with P pattern = " << highest_score << endl;                                                                                         
                        if(highest_score == highest_score_overall)
                        {
                            cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                                 << highest_score_overall << endl;
                        }                                                                             
                        if(highest_score > highest_score_overall)
                        {
                            highest_score_overall      = highest_score;
                            insertion_site_ctg_overall = insertion_site_ctg;
                            phase_assign               = false;
                        }
                        // flipped                                                          
                        long highest_score_flipped;
                        rely_on_paired = false;
                        string insertion_site_ctg_flipped = find_best_map_insertion_site_del_like(rely_on_paired,
                                                                                         this_raw_mkr,
                                                                                         this_raw_pat_flipped, 
                                                                                         this_gm_contigs,
                                                                                         contigMarkerPM,
                                                                                         &highest_score_flipped);
                        cout << "       : score with M pattern = " << highest_score_flipped << " (flipped) "<< endl;
                        if(highest_score_flipped == highest_score_overall)
                        {
                            cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                                 << highest_score_overall << endl;
                        }                                                                             
                        if(highest_score_flipped > highest_score_overall)
                        {
                            highest_score_overall      = highest_score_flipped;
                            insertion_site_ctg_overall = insertion_site_ctg_flipped;
                            phase_assign               = true;              
                        }
                        //
                        vector<string>::iterator ins_itr = std::find(contig_genetic_map[this_gmap_file].begin(),
                                                                     contig_genetic_map[this_gmap_file].end(),
                                                                     insertion_site_ctg_overall);              
                        // note del-like can be multiple for the same contig, so "sta\tend" as the label. 
                        if(phase_assign == true)
                        {                         
                            // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat_flipped)); 
                            cout << this_raw_pat_flipped     << "\t" << ss.str() << endl
                                 << contigMarkerPM[*ins_itr] << "\t" << *ins_itr << endl;
                            this_raw_pat = this_raw_pat_flipped;
                        }else
                        {
                            // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat));
                            cout << this_raw_pat             << "\t" << ss.str() << endl
                                 << contigMarkerPM[*ins_itr] << "\t" << *ins_itr << endl;                                
                        } 
                        // update genetic map and related information
                        cout << "       : inserting before " << *ins_itr
                             << " with score "               << highest_score_overall 
                             << " to  genetic map group "    << this_gmap_file
                             << " phase={"                   << phase_assign << "}"
                             << endl;
                        cout << "       : related linkage file: " << matched_map_phase[this_gmap_file] 
                             << ". " << endl << endl;
                        // target contigs in group: contig_phasing[this_linkage_file]
                        this_del_like[this_del_sta].gmap      = this_gmap_file;
                        this_del_like[this_del_sta].status    = "phased.hap.case2.2";
                        this_del_like[this_del_sta].pmpat     = this_raw_pat; // might have become "M" pattern
                        this_del_like[this_del_sta].insmarker = *ins_itr;
                        this_del_like[this_del_sta].score     = highest_score_overall; // why not highest_score_overall?
                        //
                        if(gmap_count.find(this_gmap_file) == gmap_count.end())
                        {
                            gmap_count.insert(std::pair<string, int>(this_gmap_file, 1));
                            gmap_score.insert(std::pair<string, int>(this_gmap_file, 2*this_raw_pat.size()));
                            ////gmap_start.insert(std::pair<string, unsigned long>(this_gmap_file, this_del_sta));
                        }else
                        {
                            gmap_count[this_gmap_file] += 1;
                            if(gmap_score[this_gmap_file] < 2*this_raw_pat.size())
                            {
                                gmap_score[this_gmap_file] = 2*this_raw_pat.size();
                                ////gmap_start[this_gmap_file] = this_del_sta;
                            }
                        }
                        //     
                        if(mark_score.find(this_gmap_file) == mark_score.end())
                        {
                            map<string, int> tmpscore;
                            tmpscore.insert(std::pair<string, int>(*ins_itr, highest_score_overall));
                            mark_score.insert(std::pair<string, map<string, int> >(this_gmap_file, 
                                                                                   tmpscore));
                            map<string, string> tmppmpat;
                            tmppmpat.insert(std::pair<string, string>(*ins_itr + "\t" + this_ctg, this_raw_pat));
                            mark_score_pmpat.insert(std::pair<string, map<string, string> >(this_gmap_file, 
                                                                                            tmppmpat));                                                                                        
                        }else
                        {
                            if(mark_score[this_gmap_file].find(*ins_itr) == mark_score[this_gmap_file].end())
                            {
                                mark_score[this_gmap_file].insert(std::pair<string, int>(*ins_itr, highest_score_overall));
                                mark_score_pmpat[this_gmap_file].insert(std::pair<string, string>(*ins_itr + "\t" + this_ctg, this_raw_pat));                                
                            }else
                            {
                                if(highest_score_overall > mark_score[this_gmap_file][*ins_itr])
                                {
                                    mark_score[this_gmap_file][*ins_itr] = highest_score_overall;
                                    mark_score_pmpat[this_gmap_file][*ins_itr + "\t" + this_ctg] = this_raw_pat;                                    
                                }
                            }
                        }
                        //
                        case2_hap_cnt ++;
                    }
                }
            }else
            {
                in_gmap_thus_in_phase ++;
                if(ctg_case123[3].find(this_ctg) == ctg_case123[3].end())
                {
                   ctg_case123[3].insert(std::pair<string, int>(this_ctg, 1));
                }
                // case 3: in genetic map/ordering (so in grouping/phasing) -- only phasing needed.
                string this_gmap_file = "";
                if(contig_gmap_group.find(this_left_id)  != contig_gmap_group.end())
                {
                    this_gmap_file = contig_gmap_group[this_left_id];
                }else
                if(contig_gmap_group.find(this_right_id) != contig_gmap_group.end())
                {
                    this_gmap_file = contig_gmap_group[this_right_id];
                }else;         
                vector<string> this_gm_contigs = contig_genetic_map[this_gmap_file];                                   
                cout << "       : case 3 - contig of this del-like marker \"" << ss.str() 
                     << "\" IN genetic map: "                  << this_gmap_file << endl;
                if(this_node.haphom.find("hom") != std::string::npos)
                {
                    // case 3.1: hom   
                    cout << "                  ";                                          
                    cout << "but its depth>=hom => no further genetic map completing can be done with it; "
                         << "pacbio reads overlapping such a region would go to both genotypes of group "
                         << matched_map_phase[this_gmap_file]
                         << endl;
                    this_del_like[this_del_sta].gmap   = this_gmap_file;
                    this_del_like[this_del_sta].status = "unphased.hom.case3.1"; 
                    this_del_like[this_del_sta].score  = 2*this_raw_pat.size();                                            
                    //
                    if(gmap_count.find(this_gmap_file) == gmap_count.end())
                    {
                        gmap_count.insert(std::pair<string, int>(this_gmap_file, 1));
                        gmap_score.insert(std::pair<string, int>(this_gmap_file, 2*this_raw_pat.size()));  
                        ////gmap_start.insert(std::pair<string, unsigned long>(this_gmap_file, this_del_sta));
                    }else
                    {
                        gmap_count[this_gmap_file] += 1;
                        if(gmap_score[this_gmap_file] < 2*this_raw_pat.size())
                        {
                            gmap_score[this_gmap_file] = 2*this_raw_pat.size();
                            ////gmap_start[this_gmap_file] = this_del_sta;                            
                        }
                    }
                    case3_hom_cnt ++;
                }else
                {
                    // case 3.2: hap
                    cout << "                  ";
                    cout << "and it is hap => phasing it according PM similarity "
                         << "to " << this_ctg << " in group "                
                         << matched_map_phase[this_gmap_file]
                         << " to decide genotype of pacbio reads. "
                         << endl;
                    //cout << "                  note: if " << this_ctg << " does not have enough effective pollen; "
                    //     << "using other PM of other contigs in the same group. "
                    //     << endl; // not used yet!
                    // match score to contig-itself: orignal del-like
                    assert(contigMarkerPM.find(this_left_id)  != contigMarkerPM.end());   
                    assert(contigMarkerPM.find(this_right_id) != contigMarkerPM.end());               
                    string left_ctgpat  = contigMarkerPM[this_left_id];
                    string right_ctgpat = contigMarkerPM[this_right_id];
                    //                    
                    cout << "                  Comp: investigating genetic map " << this_gmap_file << endl;
                    // find insertion site of current marker in the genetic map
                    long highest_score         = 0;
                    long highest_score_overall = 0;
                    bool phase_assign                   = false;                            
                    highest_score += find_PM_pattern_match(this_raw_pat, left_ctgpat, false);
                    highest_score += find_PM_pattern_match(this_raw_pat, right_ctgpat, false);                    
                    cout << "                       : score with P pattern to " << this_left_id << ": " << highest_score << endl;
                    if(highest_score == highest_score_overall)
                    {
                        cout << "                   Warning: a tie found for " << this_raw_mkr << " with score " 
                             << highest_score_overall << endl;
                    }                                                                             
                    if(highest_score > highest_score_overall)
                    {
                        highest_score_overall      = highest_score;
                        phase_assign               = false;
                    }
                    // flipped                                                          
                    long highest_score_flipped = 0;
                    highest_score_flipped += find_PM_pattern_match(this_raw_pat_flipped, left_ctgpat, false);
                    highest_score_flipped += find_PM_pattern_match(this_raw_pat_flipped, right_ctgpat, false);                    
                    cout << "                       : score with M pattern to " << this_left_id << ": " << highest_score_flipped
                         << " (flipped) " << endl;                                        
                    if(highest_score_flipped == highest_score_overall)
                    {
                        cout << "                   Warning: a tie found for " << this_raw_mkr << " with score " 
                             << highest_score_overall << endl;
                    }                                                                             
                    if(highest_score_flipped > highest_score_overall)
                    {
                        highest_score_overall      = highest_score_flipped;
                        phase_assign               = true;              
                    }          
                    // note del-like can be multiple for the same contig, so "sta\tend" as the label. 
                    if(phase_assign == true)
                    {                         
                        // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat_flipped)); 
                        cout << this_raw_pat_flipped         << "\t" << ss.str()     << endl
                             << contigMarkerPM[this_left_id] << "\t" << this_left_id << endl;
                        this_raw_pat = this_raw_pat_flipped;
                    }else
                    {
                        // contigMarkerPM.insert(std::pair<string, string>(ss.str(), this_raw_pat));
                        cout << this_raw_pat                 << "\t" << ss.str()     << endl
                             << contigMarkerPM[this_left_id] << "\t" << this_left_id << endl;
                    }
                    // update genetic map and related information
                    cout << "                  : no need to insert - " << this_left_id
                         << " with score "                  << highest_score_overall 
                         << " to  genetic map group "       << this_gmap_file
                         << " phase={"                      << phase_assign << "}"
                         << endl;
                    cout << "                  : related linkage file: " << matched_map_phase[this_gmap_file] 
                         << ". " << endl << endl;
                    //
                    this_del_like[this_del_sta].gmap      = this_gmap_file; // need to update with the insertion site
                    this_del_like[this_del_sta].status    = "phased.hap.case3.2";      
                    this_del_like[this_del_sta].pmpat     = this_raw_pat;   // might have become "M" pattern         
                    this_del_like[this_del_sta].insmarker = this_right_id; 
                    this_del_like[this_del_sta].score     = highest_score_overall;                                                                  
                    //
                    if(gmap_count.find(this_gmap_file) == gmap_count.end())
                    {
                        gmap_count.insert(std::pair<string, int>(this_gmap_file, 1));
                        gmap_score.insert(std::pair<string, int>(this_gmap_file, 2*this_raw_pat.size()));
                        ////gmap_start.insert(std::pair<string, unsigned long>(this_gmap_file, this_del_sta));
                    }else
                    {
                        gmap_count[this_gmap_file] += 1;
                        if(gmap_score[this_gmap_file] < 2*this_raw_pat.size())
                        {
                            gmap_score[this_gmap_file] = 2*this_raw_pat.size();
                            ////gmap_start[this_gmap_file] = this_del_sta;                            
                        }                        
                    } 
                    case3_hap_cnt ++;                                                                                                                                                      
                }
            }
            // next contig-region.
            ditr ++;
        }
        // find out if we can assign a del-like marker to a genetic map / phasing group 
        // same contig-del-like should all go to the same group:
        cout << "         Genetic map insertion cases summary for " << this_ctg << " related del-like: " << endl;
        cout << "            checked " << this_del_like.size() << " del-like regions: " << endl;
        map<string, int>::iterator gmcitr;
        map<string, int>::iterator gmcitr_end;
        gmcitr     = gmap_count.begin();
        gmcitr_end = gmap_count.end();
        int    max          = 0;
        string max_gmapfile = "gmap.file.unavailable";
        int    max_cnt      = 0;
        while(gmcitr != gmcitr_end)
        {
            // do not score more than o for unvailable gmap
            if((*gmcitr).first.find("unavailable") != std::string::npos && max_cnt > 0)
            {
                gmcitr ++;
                continue;
            }
            //
            cout << "               inserted into genetic map ";
            cout << (*gmcitr).first             << "\t" 
                 << (*gmcitr).second            << " times with max-score " 
                 << gmap_score[(*gmcitr).first] << endl;
            if((*gmcitr).second > max)
            {
                max_gmapfile = (*gmcitr).first;
                max          = (*gmcitr).second;
                max_cnt      = 1;
            }else 
            if((*gmcitr).second == max)
            {
                if(gmap_score[max_gmapfile] < gmap_score[(*gmcitr).first])
                {                 
                    cout << "            ";
                    cout << "equal occurrence but grouping max-score " << gmap_score[(*gmcitr).first] 
                         << " with " << (*gmcitr).first
                         << " is higher than last score "              << gmap_score[max_gmapfile]    
                         << " with " << max_gmapfile 
                         << ", so max_gmapfile will be updated as " 
                         << (*gmcitr).first << endl;
                    max_gmapfile = (*gmcitr).first;
                    max          = (*gmcitr).second;                           
                }else
                if(gmap_score[max_gmapfile] == gmap_score[(*gmcitr).first])
                {
                    max_cnt ++;                    
                }else ;
            }else
            {
                if(max_gmapfile.compare("gmap.file.unavailable") == 0 && 
                   (*gmcitr).first.compare("gmap.file.unavailable") != 0)
                {
                    max_gmapfile = (*gmcitr).first;
                    max          = (*gmcitr).second;
                    max_cnt      = 1;                    
                }
            }                      
            gmcitr ++;
        }
        if(max_cnt > 1)
        {
            cout << "         current max observed more than once with same score. " << endl;
        }
        cout << "         this contig would go to genetic map "                      << max_gmapfile 
             << " (if it has not been collected in current genetic map formed by snp-contig-markers)." << endl;
        // find out, for max_gmapfile, which snp-contig-marker should be inserted.        
        map<string, map<string, int> >::iterator gmapitr = mark_score.find(max_gmapfile); // <gmapfile, <contig_marker, insertion_score> >
        if(gmapitr == mark_score.end())
        {
                cout << "            ";
                cout << "This del-like cannot be inserted into genetic map (or it is already in genetic map). "
                     << endl;
        }else
        {
            map<string, int> this_mark_score = (*gmapitr).second;
            map<string, int>::iterator scoitr;
            map<string, int>::iterator scoitr_end;
            scoitr     = this_mark_score.begin();
            scoitr_end = this_mark_score.end();           
            int    max_score  = 0;
            string max_marker = ""; 
            int    max_score_cnt = 0;
            string max_pmpat  = "";
            while(scoitr != scoitr_end)
            {
                if( (*scoitr).second > max_score)
                {
                    max_score  = (*scoitr).second;
                    max_marker = (*scoitr).first;
                    max_score_cnt = 1;
                    max_pmpat  = mark_score_pmpat[max_gmapfile][max_marker + "\t" + this_ctg];
                }else
                if( (*scoitr).second == max_score)                
                {
                    max_score_cnt ++;
                }else ;
                scoitr ++;
            }
            assert(mark_score.find(max_gmapfile) != mark_score.end());
            assert(mark_score[max_gmapfile].find(max_marker) != mark_score[max_gmapfile].end());                
            cout << "            ";
            cout << "Insert this contig " << this_ctg << " before " << max_marker
                 << " of " << max_gmapfile << " with original score " 
                 << mark_score[max_gmapfile][max_marker] << endl;
            assert(to_be_inserted_del_contig_id.find(this_ctg) == to_be_inserted_del_contig_id.end());     
            to_be_inserted_del_contig_id.insert(std::pair<string, string>(this_ctg, max_gmapfile+"#"+max_marker));
            to_be_inserted_del_contig_pmpat.insert(std::pair<string, string>(this_ctg, max_pmpat));     
             
            // update genetic map - 2020-01-25
            
            /*
            cout  << "   check: updating genetic map " << max_gmapfile << " with " << this_left_id << endl;
            
            assert((*contig_genetic_map_updated).find(max_gmapfile) != (*contig_genetic_map_updated).end());            
            vector<string>::iterator tmp_ins_itr = std::find((*contig_genetic_map_updated)[max_gmapfile].begin(),
                                                             (*contig_genetic_map_updated)[max_gmapfile].end(),
                                                              max_marker);  
            cout << "   check: insert before " << *tmp_ins_itr << " done. " << endl;
            (*contig_genetic_map_updated)[max_gmapfile].insert(tmp_ins_itr, this_left_id+"_del");
            (*contig_genetic_map_updated)[max_gmapfile].insert(tmp_ins_itr, this_right_id+"_del");   
            cout << "   check: contig_genetic_map_updated. " << endl;         
            // update genetic map
            cout << "   check: number of markers in (*contig_gmap_group_updated) = " << (*contig_gmap_group_updated).size() << endl;
            (*contig_gmap_group_updated).insert(std::pair<string, string>(this_left_id+"_del",  max_gmapfile));              
            (*contig_gmap_group_updated).insert(std::pair<string, string>(this_right_id+"_del", max_gmapfile));
            cout << "   check: contig_gmap_group_updated. " << endl;                     
            // update contig-marker-PM map
            contigMarkerPM_updated.insert(std::pair<string, string>(this_left_id+"_del",  max_pmpat)); 
            contigMarkerPM_updated.insert(std::pair<string, string>(this_right_id+"_del", max_pmpat));                
            cout << "   check: updating genetic map " << max_gmapfile  << " done. " << endl;            
            */
        }
        // update this_del_like again according to final decision; update phasing status
        cout << "            Some more details about insertion of the del-like marker: " << endl;                
        map<unsigned long, DELMARKER>::iterator tdlitr;
        map<unsigned long, DELMARKER>::iterator tdlitr_end;
        tdlitr     = this_del_like.begin();
        tdlitr_end = this_del_like.end();
        while(tdlitr != tdlitr_end)
        {
            string this_phase = "P";
            string this_linkage_file = "unphased.txt";
            if(matched_map_phase.find(max_gmapfile) != matched_map_phase.end())
            {
                this_linkage_file = matched_map_phase[max_gmapfile];
            }
            if((*tdlitr).second.pmpat.find("P")==std::string::npos) this_phase = "M";            
            ofp  << (*tdlitr).second.ctg            << "\t" 
                 << (*tdlitr).second.sta            << "\t" 
                 << (*tdlitr).second.end            << "\t" 
                 << (*tdlitr).second.haphom         << "\t" 
                 << this_phase                      << "\t" 
                 << (*tdlitr).second.score          << "\t" 
                 << this_linkage_file               << "\t"
                 << max_gmapfile                    << "\t"
                 << (*tdlitr).second.pmpat          << endl;
            std::stringstream ssdel;
            ssdel.str("");
            ssdel << (*tdlitr).second.ctg << ":" << (*tdlitr).second.sta << "-" << (*tdlitr).second.end;
            // get position to insert in genetic map
            if((*tdlitr).second.status.find("unphased.hom.case1.1")   != std::string::npos)
            {
                cout << "               ";                
                cout << "This del-like " << ssdel.str() << " in hom regions - not phased - go to unmapped. " << endl;
            }else
            if((*tdlitr).second.status.find("unphased.hap.case1.2.1") != std::string::npos)
            {
                cout << "               ";                
                cout << "This del-like " << ssdel.str() << " in hap regions - not phased - go to unmapped. " << endl;
            }else  
            if((*tdlitr).second.status.find("phased.hap.case1.2.2")   != std::string::npos)
            {
                assert(mark_score.find((*tdlitr).second.gmap) != mark_score.end());
                assert(mark_score[(*tdlitr).second.gmap].find((*tdlitr).second.insmarker) !=
                       mark_score[(*tdlitr).second.gmap].end());
                cout << "               ";                                
                cout << "Insert this del-like " << ssdel.str() << " before " << (*tdlitr).second.insmarker 
                     << " of " << max_gmapfile << " with original score " 
                     << mark_score[(*tdlitr).second.gmap][(*tdlitr).second.insmarker] << endl;                
            }else
            if((*tdlitr).second.status.find("unphased.hom.case2.1")   != std::string::npos)
            {
                cout << "               ";                
                cout << "This del-like " << ssdel.str() << " in hom regions - not phased (but snp-marker phased). "<< endl;
            }else                        
            if((*tdlitr).second.status.find("phased.hap.case2.2")     != std::string::npos)            
            {
                assert(mark_score.find((*tdlitr).second.gmap) != mark_score.end());
                assert(mark_score[(*tdlitr).second.gmap].find((*tdlitr).second.insmarker) !=
                       mark_score[(*tdlitr).second.gmap].end());                
                cout << "               ";                                
                cout << "Insert this del-like " << ssdel.str() << " before " << (*tdlitr).second.insmarker 
                     << " of " << max_gmapfile << " with original score " 
                     << mark_score[(*tdlitr).second.gmap][(*tdlitr).second.insmarker] << endl;
            }else
            if((*tdlitr).second.status.find("unphased.hom.case3.1")   != std::string::npos)
            {
                cout << "               ";
                cout << "This del-like " << ssdel.str() << " in hom regions - not phased (but in genetic map). " << endl;
            }else            
            if((*tdlitr).second.status.find("phased.hap.case3.2")     != std::string::npos)
            {
                cout << "               ";
                cout << "This del-like " << ssdel.str() << " in hap regions - and phased (but in genetic map). " << endl;
            } 
            else ;            
            //
            tdlitr ++;
        }
        cout << endl;
        //
        (*delctgitr).second.clear();
        (*delctgitr).second.insert(this_del_like.begin(), this_del_like.end());
        // next contig
        delctgitr ++;
    }  
    // del-like region summary      
    // all phased are in genetic map
    //  because phased ones not in genetic map have been inserted into genetic maps by comparing snp-contig PM patterns.
    cout << "   Info: analyzed " << caseall_total << " del-like markers, including: "              << endl
         << "         case1 ......... no  genetic map + no  phasing: " << not_in_gmap_not_in_phase << " from "   
         << ctg_case123[1].size() << " contigs. "                                                  << endl
         << "         case1_bad_cnt ...............................: " << case1_bad_cnt            << endl
         << "         case1_hap_cnt ...............................: " << case1_hap_cnt            << endl
         << "         case1_hom_cnt ...............................: " << case1_hom_cnt            << endl
         << "         case2 ......... no  genetic map + yes phasing: " << not_in_gmap_but_in_phase << " from "  
         << ctg_case123[2].size() << " contigs. "                                                  << endl              
         << "         case2_hap_cnt ...............................: " << case2_hap_cnt            << endl
         << "         case2_hom_cnt ...............................: " << case2_hom_cnt            << endl
         << "         case3 ..........yes genetic map + yes phasing: " << in_gmap_thus_in_phase    << " from "
         << ctg_case123[3].size() << " contigs. "                                                  << endl
         << "         case3_hap_cnt ...............................: " << case3_hap_cnt            << endl
         << "         case3_hom_cnt ...............................: " << case3_hom_cnt            << endl;
    //
    ofp.close(); // phased file
    ///////////////////////////////////////// finalizing ///////////////////////////////////////////////////////////////
    //////////////// update genetic map with del-contig markers (only those not in existing genetic map) ///////////////
    map<string, vector<string> > tmp_gm;    
    tmp_gm.insert(contig_genetic_map.begin(), contig_genetic_map.end());
    map<string, string>          tmp_group;
    tmp_group.insert(contig_gmap_group.begin(), contig_gmap_group.end());    
    map<string, string>          tmp_pmpat;
    tmp_pmpat.insert(contigMarkerPM.begin(), contigMarkerPM.end());
    //
    map<string, string>::iterator tbitr;
    map<string, string>::iterator tbitr_end;
    tbitr     = to_be_inserted_del_contig_id.begin();
    tbitr_end = to_be_inserted_del_contig_id.end();
    cout << endl;
    cout << "   Info: " << to_be_inserted_del_contig_id.size() 
         << " del-contigs need to be inserted in genetic map. " << endl;
    while(tbitr != tbitr_end)
    {    
        string this_ctg      = (*tbitr).first;      // del-related contig not in genetic map yet
        string this_left_id  = this_ctg + "\tleft"; // building ids like ctg\tleft/right in asPollinator format
        string this_right_id = this_ctg + "\tright";
        vector<string> gm_mkr_info = split_string((*tbitr).second, '#');   
        assert(gm_mkr_info.size()==2);         
        string max_gmapfile = gm_mkr_info[0];    
        string max_marker   = gm_mkr_info[1];     
        map<string, vector<string> >::iterator gmvitr = tmp_gm.find(max_gmapfile);           
        assert( gmvitr!= tmp_gm.end());                    
        vector<string>::iterator tmp_ins_itr = std::find( ((*gmvitr).second).begin(),
                                                          ((*gmvitr).second).end(),
                                                          max_marker);  
        assert(tmp_ins_itr != ((*gmvitr).second).end());                                                                                                            
        // this results in malloc error?
        ////vector<string>::iterator itr1 = ((*gmvitr).second).insert(tmp_ins_itr, this_left_id+"_del");  
        ////vector<string>::iterator itr2 = ((*gmvitr).second).insert(tmp_ins_itr, this_right_id+"_del");           
        vector<string> tmpgmcontig;
        vector<string>::iterator tmpcitr;
        vector<string>::iterator tmpcitr_end;
        tmpcitr     = (*gmvitr).second.begin();
        tmpcitr_end = (*gmvitr).second.end();
        while(tmpcitr != tmpcitr_end)
        {
            if((*tmpcitr).compare(max_marker) == 0)
            {
                tmpgmcontig.push_back(this_left_id+"_del");
                tmpgmcontig.push_back(this_right_id+"_del");
            }
            tmpgmcontig.push_back(*tmpcitr);
            tmpcitr ++;
        }
        (*gmvitr).second.clear();
        (*gmvitr).second = tmpgmcontig;
        std::copy ( tmpgmcontig.begin(), tmpgmcontig.end(), (*gmvitr).second.begin() );        
        // update genetic map
        tmp_group.insert(std::pair<string, string>(this_left_id+"_del",  max_gmapfile));              
        tmp_group.insert(std::pair<string, string>(this_right_id+"_del", max_gmapfile));        
        // update contig-marker-PM map
        string max_pmpat = to_be_inserted_del_contig_pmpat[this_ctg];
        tmp_pmpat.insert(std::pair<string, string>(this_left_id+"_del",  max_pmpat)); 
        tmp_pmpat.insert(std::pair<string, string>(this_right_id+"_del", max_pmpat));                        
        // 
        tbitr ++;
    }  
    // update final genetic map etc
    (*contig_genetic_map_updated).insert(tmp_gm.begin(), tmp_gm.end());      // final genetic map
    (*contig_gmap_group_updated).insert(tmp_group.begin(), tmp_group.end()); // not used further
    contigMarkerPM_updated.insert(tmp_pmpat.begin(), tmp_pmpat.end());       // final pm pattern (with flipping)
    // output final genetic map corrected with snp and del markers for later scaffolding
    if(1)
    {
        // create an intermediate folder for collecting updated genetic maps
        std::stringstream iss;
        iss.str("");
        iss << "./z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs_" << hap_low_to_hom_low_reset;
        string tmpfolder = iss.str();
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
        unsigned long total_gm_contig_markers = 0;
        map<string, int> total_gm_contig_uniq;
        map<string, vector<string> >::iterator upd_mapfitr;
        map<string, vector<string> >::iterator upd_mapfitr_end;
        upd_mapfitr     = (*contig_genetic_map_updated).begin();
        upd_mapfitr_end = (*contig_genetic_map_updated).end();
        while(upd_mapfitr != upd_mapfitr_end)
        {
            string         this_gm_file    = (*upd_mapfitr).first;
            vector<string> this_gm_contigs = (*upd_mapfitr).second;
            cout << "   Check: map "       << this_gm_file       << " update to " 
                 << this_gm_contigs.size() << " contig markers: "<< endl;                 
            total_gm_contig_markers += this_gm_contigs.size();
            //
            string         this_gm_file_upd   = tmpfolder + "/upd_"   + this_gm_file; // genetic map
            string         this_gm_file_updPM = tmpfolder + "/PM_upd" + this_gm_file; // PM values of genetic map
            ofstream updofp;
            updofp.open(this_gm_file_upd.c_str(), ios::out);
            if(!updofp.good())
            {
                cout << "   Error: cannot open file to collected updated genetic map: " 
                     << this_gm_file_upd   << endl;
                return false;
            }
            ofstream updPMofp;
            updPMofp.open(this_gm_file_updPM.c_str(), ios::out);
            if(!updPMofp.good())
            {
                cout << "   Error: cannot open file to collected updated PMs of genetic map: " 
                     << this_gm_file_updPM << endl;
                return false;
            }
            //       
            vector<string>::iterator gmcitr     = this_gm_contigs.begin();
            vector<string>::iterator gmcitr_end = this_gm_contigs.end();
            while(gmcitr != gmcitr_end)
            {
                // 
                vector<string> markerinfo = split_string(*gmcitr, '\t'); // contig\tleft/right or contig\tcontig
                if((*gmcitr).find("left") == std::string::npos && 
                   (*gmcitr).find("right")==std::string::npos )
                {
                    if(total_gm_contig_uniq.find(markerinfo[1]) == total_gm_contig_uniq.end())
                    {
                        total_gm_contig_uniq.insert(std::pair<string, int>(markerinfo[1], 1));
                    }
                }else
                {
                    if(total_gm_contig_uniq.find(markerinfo[0]) == total_gm_contig_uniq.end())
                    {
                        total_gm_contig_uniq.insert(std::pair<string, int>(markerinfo[0], 1));
                    }                    
                }
                // genetic map
                cout   << "        " << *gmcitr << endl;
                updofp << *gmcitr;
                map<string, GNMAP>::iterator disitr = gmContig.find(*gmcitr);
                if(disitr != gmContig.end())
                {
                    updofp << "\t" << (*disitr).second.distance;
                }else
                {
                    updofp << "\t" << "NA";
                }
                if(gmcitr!= gmcitr_end-1)
                {
                    double new_rcomb_freq = calculate_recomb_freq(contigMarkerPM_updated[*gmcitr],
                                                                  contigMarkerPM_updated[*(gmcitr+1)]);
                    updofp << "\t" << new_rcomb_freq;
                }
                else
                {
                    updofp << "\t" << "NA";
                }
                updofp << endl;
                // PM
                updPMofp << contigMarkerPM_updated[*gmcitr] << "\t" << *gmcitr << endl;               
                //
                gmcitr ++;
            }
            //
            updofp.close();
            updPMofp.close();
            // next map
            upd_mapfitr ++;
        }
        cout << "   Info: final number of contig markers in genetic map: "  << total_gm_contig_markers << endl
             << "         corresponding to " << total_gm_contig_uniq.size() << " unique contigs."      << endl;
    }       
    //
    return true;
}
// complete 2: complete genetic map with snp-contig-markers
bool complete_genetic_map_snp_based(map<string, string>             contigMarkerPM, 
                                    map<string, map<string, bool> > contig_phasing,
                                    map<string, string>             ctgGroup,
                                    map<string, vector<string> >    contig_genetic_map, 
                                    map<string, string>             contig_gmap_group,                        
                                    map<string, string>             matched_phase_map,
                                    map<string, string>             matched_map_phase,
                                    map<string, GNMAP>              gmContig,
                                    map<string, vector<string> >*   contig_genetic_map_updated,
                                    map<string, bool>*              ctgJoinMapPhasing_all,  
                                    map<string, string>*            ctgGroup_updated,
                                    map<string, string>*            contig_gmap_group_updated,
                                    map<string, string>*            contigMarkerPM_updated
                                    )
{
    /*
    allraw contig-left/right markers           - contigMarkerPM       : <contig\tleft/right, marker-PM-pattern>
    -- grouping related
    phased contig-left/right markers in groups - contig_phasing       : <this_linkage_file , <contig\tleft/right, bool_flip> > 
    phased contig-left/right markers           - ctgGroup             : <contig\tleft/right, this_linksage_file>
    contig phasing information                 - ctgJoinMapPhasing_all: <contig\tleft/right, bool_flip>    
    -- ordering related
    mapped contig-left/right markers           - contig_genetic_map   : <this_gmap_file    , <contig\tleft/right> >
    mapped contig-left/right markers in maps   - contig_gmap_group    : <contig\tleft/right, this_gmap_file> 
    -- 
    matching info of phased file and map file  - matched_phase_map    : <this_linkage_file , this_gmap_file>  
    */
    // initialize the genetic map with mapped contig markers
    (*contig_genetic_map_updated).clear();
    *ctgGroup_updated          = ctgGroup;
    *contig_gmap_group_updated = contig_gmap_group;
    *contigMarkerPM_updated    = contigMarkerPM;
    // ctgJoinMapPhasing_all: this includes initial phased contigs; to be updated by integrating ungrouped_ctgPhasing
    map<string, bool> ungrouped_ctgPhasing;
    //
    map<string, string>::iterator rawmkritr;
    map<string, string>::iterator rawmkritr_end;
    rawmkritr     = contigMarkerPM.begin();
    rawmkritr_end = contigMarkerPM.end();
    while(rawmkritr != rawmkritr_end)
    {
        string this_raw_mkr = (*rawmkritr).first;  // ctg\tleft/right
        string this_raw_pat = (*rawmkritr).second; // PM pattern: already flipped if phase info from JoinMap agreed.
        string this_raw_pat_flipped = flip_PMpat(this_raw_pat);        
        // level 1: check if current marker in genetic map or not = ordering checking           
        if(contig_gmap_group.find(this_raw_mkr) == contig_gmap_group.end())
        {
            cout << "   Info: this marker " << this_raw_mkr << " not in genetic map: " << endl;
            // level 2: check if current marker in phased group or not = phasing checking
            if(ctgGroup.find(this_raw_mkr) == ctgGroup.end())
            {
                // check current contig marker with  all phased groups
                // all patterns not in contigMarkerPM has not been corrected with phasing 
                //     -- need check both initial and flipped and reverse.
                cout << "       : this marker "                 << this_raw_mkr 
                     << " not in phased groups." << endl;
                //
                map<string, vector<string> >::iterator gmfitr;
                map<string, vector<string> >::iterator gmfitr_end;
                gmfitr     = contig_genetic_map.begin();
                gmfitr_end = contig_genetic_map.end();
                long highest_score_overall = 0;
                string insertion_site_ctg_overall   = "";
                string this_gmap_file_overall       = "";
                bool phase_assign                   = false; 
                while(gmfitr != gmfitr_end)
                {
                    string         this_gmap_file   = (*gmfitr).first;
                    vector<string> this_gm_contigs  = (*gmfitr).second;
                    // raw
                    long highest_score;
                    bool rely_on_paired = false;
                    string insertion_site_ctg = find_best_map_insertion_site(rely_on_paired,
                                                                             this_raw_mkr,
                                                                             this_raw_pat, 
                                                                             this_gm_contigs,
                                                                             contigMarkerPM,
                                                                             &highest_score);
                    if(highest_score == highest_score_overall)
                    {
                        cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                             << highest_score_overall << endl;
                    }                                                                             
                    if(highest_score > highest_score_overall)
                    {
                        highest_score_overall      = highest_score;
                        insertion_site_ctg_overall = insertion_site_ctg;
                        this_gmap_file_overall     = this_gmap_file;
                        phase_assign               = false;
                    }
                    // flipped
                    long highest_score_flipped;
                    rely_on_paired = false;                    
                    string insertion_site_ctg_flipped = find_best_map_insertion_site(rely_on_paired,
                                                                             this_raw_mkr,
                                                                             this_raw_pat_flipped, 
                                                                             this_gm_contigs,
                                                                             contigMarkerPM,
                                                                             &highest_score_flipped);
                    if(highest_score_flipped == highest_score_overall)
                    {
                        cout << "         Warning: a tie found for " << this_raw_mkr << " with score " 
                             << highest_score_overall << endl;
                    }                                                                             
                    if(highest_score_flipped > highest_score_overall)
                    {
                        highest_score_overall      = highest_score_flipped;
                        insertion_site_ctg_overall = insertion_site_ctg_flipped;
                        this_gmap_file_overall     = this_gmap_file;          
                        phase_assign               = true;              
                    }                                                          
                    //
                    gmfitr ++;
                }
                vector<string>::iterator ins_itr = std::find(contig_genetic_map[this_gmap_file_overall].begin(),
                                                             contig_genetic_map[this_gmap_file_overall].end(),
                                                             insertion_site_ctg_overall);
                // update genetic map and related information
                cout << "       : inserting before " << *ins_itr
                     << " with score "               << highest_score_overall 
                     << " to  genetic map group "    << this_gmap_file_overall 
                     << " phase={"                   << phase_assign << "}"
                     << endl;                
                contig_genetic_map[this_gmap_file_overall].insert(ins_itr, this_raw_mkr);   
                // new phase status for ungrouped contigs
                ungrouped_ctgPhasing.insert(std::pair<string, bool>(this_raw_mkr, phase_assign));    
                // update
                (*contig_gmap_group_updated).insert(std::pair<string, string>(this_raw_mkr, this_gmap_file_overall));
                // update
                assert(matched_map_phase.find(this_gmap_file_overall) != matched_map_phase.end());
                cout << "       : related PM phasing file: " << matched_map_phase[this_gmap_file_overall] << ". " << endl;                 
                (*ctgGroup_updated).insert(std::pair<string, string>(this_raw_mkr, 
                                                                     matched_map_phase[this_gmap_file_overall] ));
                //update
                if(phase_assign == true)
                {                    
                    (*contigMarkerPM_updated).insert(std::pair<string, string>(this_raw_mkr, this_raw_pat_flipped));
                }else
                {
                    (*contigMarkerPM_updated).insert(std::pair<string, string>(this_raw_mkr, this_raw_pat));                    
                }
                cout << "   Info: new PM pattern collected for " << this_raw_mkr << endl;
            }else
            {
                // check current contig marker within one phased group                
                // all patterns in contigMarkerPM has been corrected with phasing, if they are in phased groups.
                //     so no flip needed for such PM patterns.
                string this_linkage_file       = ctgGroup[this_raw_mkr];
                string this_gmap_file          = matched_phase_map[this_linkage_file];
                vector<string> this_gm_contigs = contig_genetic_map[this_gmap_file];                
                cout << "       : this marker "            << this_raw_mkr 
                     << " --- in phased group "            << this_linkage_file 
                     << "; matching genetic map file "     << this_gmap_file
                     << " with " << this_gm_contigs.size() << " ctg-markers. " 
                     << endl;
                long highest_score;  
                bool rely_on_paired = true;                  
                string insertion_site_ctg = find_best_map_insertion_site(rely_on_paired,
                                                                         this_raw_mkr,
                                                                         this_raw_pat, 
                                                                         this_gm_contigs,
                                                                         contigMarkerPM,
                                                                         &highest_score);
                vector<string>::iterator ins_itr = std::find(contig_genetic_map[this_gmap_file].begin(),
                                                             contig_genetic_map[this_gmap_file].end(),
                                                             insertion_site_ctg);
                // update genetic map and related information
                cout << "       : inserting before " << *ins_itr
                     << " with score "               << highest_score 
                     << " to  genetic map group "    << this_gmap_file 
                     << endl;
                contig_genetic_map[this_gmap_file].insert(ins_itr, this_raw_mkr);   
                // update
                (*contig_gmap_group_updated).insert(std::pair<string, string>(this_raw_mkr, this_gmap_file));  
                // update: if it is in phased group, it is collected already.
                assert( (*contigMarkerPM_updated).find(this_raw_mkr) != (*contigMarkerPM_updated).end() );
            }
        }else
        {
            // cout << "   Check: " << this_raw_mkr << " in map " << contig_gmap_group[this_raw_mkr] << endl;
        }
        // next
        rawmkritr ++;
    }   
    //    
    (*contig_genetic_map_updated).insert(contig_genetic_map.begin(), contig_genetic_map.end());
    //
    if(1)
    {
        // create an intermediate folder for collecting updated genetic maps
        std::stringstream iss;
        iss.str("");
        iss << "./z_genetic_maps_updated_with_PMsimilarity_of_snp_contigs_intermediate_" << hap_low_to_hom_low_reset;
        string tmpfolder = iss.str();
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
        unsigned long total_gm_contig_markers = 0;
        map<string, vector<string> >::iterator upd_mapfitr;
        map<string, vector<string> >::iterator upd_mapfitr_end;
        upd_mapfitr     = (*contig_genetic_map_updated).begin();
        upd_mapfitr_end = (*contig_genetic_map_updated).end();
        while(upd_mapfitr != upd_mapfitr_end)
        {
            string         this_gm_file    = (*upd_mapfitr).first;
            vector<string> this_gm_contigs = (*upd_mapfitr).second;
            cout << "   Check: map "       << this_gm_file       << " update to " 
                 << this_gm_contigs.size() << " contig markers: "<< endl;                 
            total_gm_contig_markers += this_gm_contigs.size();
            //
            string         this_gm_file_upd  = tmpfolder + "/upd_"   + this_gm_file; // genetic map
            string         this_gm_file_updPM = tmpfolder + "/PM_upd" + this_gm_file; // PM values of genetic map
            ofstream updofp;
            updofp.open(this_gm_file_upd.c_str(), ios::out);
            if(!updofp.good())
            {
                cout << "   Error: cannot open file to collected updated genetic map: " 
                     << this_gm_file_upd   << endl;
                return false;
            }
            ofstream updPMofp;
            updPMofp.open(this_gm_file_updPM.c_str(), ios::out);
            if(!updPMofp.good())
            {
                cout << "   Error: cannot open file to collected updated PMs of genetic map: " 
                     << this_gm_file_updPM << endl;
                return false;
            }
            //       
            vector<string>::iterator gmcitr     = this_gm_contigs.begin();
            vector<string>::iterator gmcitr_end = this_gm_contigs.end();
            while(gmcitr != gmcitr_end)
            {
                // genetic map
                cout   << "        " << *gmcitr << endl;
                updofp << *gmcitr;
                map<string, GNMAP>::iterator disitr = gmContig.find(*gmcitr);
                if(disitr != gmContig.end())
                {
                    updofp << "\t" << (*disitr).second.distance;
                }else
                {
                    updofp << "\t" << "NA";
                }
                if(gmcitr!= gmcitr_end-1)
                {
                    double new_rcomb_freq = calculate_recomb_freq(contigMarkerPM[*gmcitr],
                                                                  contigMarkerPM[*(gmcitr+1)]);
                    updofp << "\t" << new_rcomb_freq;
                }else
                {
                    updofp << "\t" << "NA";
                }
                updofp << endl;
                // PM
                updPMofp << contigMarkerPM[*gmcitr] << "\t" << *gmcitr << endl;               
                //
                gmcitr ++;
            }
            //
            updofp.close();
            updPMofp.close();
            // next map
            upd_mapfitr ++;
        }
        cout << "   Info: final number of contig markers in genetic map: " << total_gm_contig_markers << endl;
    }
    // update the groups with additional contigs not in the initial linkage groups from JoinMap
    (*ctgJoinMapPhasing_all).insert(ungrouped_ctgPhasing.begin(), ungrouped_ctgPhasing.end());
    cout << "   Info: final phasing status collected for " << (*ctgJoinMapPhasing_all).size() << " markers. " << endl;
    //
    return true;
}
//
double calculate_recomb_freq(string PMpat1, string PMpat2)
{
    double recomb_freq = 0.0;
    string pm1 = PMpat1;
    string pm2 = PMpat2;
    std::transform(pm1.begin(), pm1.end(), pm1.begin(), ::toupper);    
    std::transform(pm2.begin(), pm2.end(), pm2.begin(), ::toupper);      
    // only consider effective exact P/M positions
    unsigned long matchedp   = 0;
    unsigned long effectivep = 0;
    assert(pm1.size() == pm2.size());
    for(unsigned long ii = 0; ii < pm1.size(); ii ++)
    {
        string letter1 = pm1.substr(ii, 1);
        string letter2 = pm2.substr(ii, 1);
        //
        if(letter1.compare("U")!=0 && letter2.compare("U") != 0)
        {
            effectivep ++;
        } 
        //
        if(letter1.compare("P")==0 && letter2.compare("P") == 0)
        {
            matchedp ++;
        }else
        if(letter1.compare("M")==0 && letter2.compare("M") == 0)
        {
            matchedp ++;
        }else ;
    }            
    recomb_freq = (effectivep-matchedp)*1.0/effectivep*100;
    //
    cout << PMpat1 << endl
         << PMpat2 << endl;
    //
    cout << "   Info: effectivep.=" << effectivep              << endl
         << "         recombined.=" << effectivep - matchedp   << endl
         << "         recomb_freq=" << recomb_freq             << endl;
    
    return recomb_freq;
}
//
string find_best_map_insertion_site_del_like(bool                rely_on_paired,
                                             string              this_raw_mkr,
                                             string              this_raw_pat, 
                                             vector<string>      this_gm_contigs,
                                             map<string, string> contigMarkerPM,
                                             long*               highest_score)
{
    // position in a genetic map before which this raw marker will be inserted
    vector<string>::iterator highest_itr;
    // checking: other markers
    vector<string>::iterator gmcitr;
    vector<string>::iterator gmcitr_end;
    gmcitr                       =    this_gm_contigs.begin();
    gmcitr_end                   = -- this_gm_contigs.end(); // last second 
    *highest_score  = 0;
    long second_highest = 0;
    highest_itr                  = gmcitr; // will insert between gmcitr and gmcitr+1
    while(gmcitr != gmcitr_end)
    {
        string left_ctgmkr  = *gmcitr;
        string right_ctgmkr = *(gmcitr+1);  
        //
        vector<string> leftinfo  = split_string(left_ctgmkr, '\t');
        vector<string> rightinfo = split_string(right_ctgmkr, '\t');
        if(leftinfo[0].compare(rightinfo[0])==0)
        {
            // do not insert between the same contig; namely 
            //    not between a\tright and b\tleft when a==b
            gmcitr ++;
            continue;
        }
        //
        assert(contigMarkerPM.find(left_ctgmkr)  != contigMarkerPM.end());   
        assert(contigMarkerPM.find(right_ctgmkr) != contigMarkerPM.end());               
        string left_ctgpat  = contigMarkerPM[left_ctgmkr];
        string right_ctgpat = contigMarkerPM[right_ctgmkr];
        //
        long matchedp = 0;
        matchedp += find_PM_pattern_match(this_raw_pat, left_ctgpat, false);
        matchedp += find_PM_pattern_match(this_raw_pat, right_ctgpat, false);
        if(matchedp > *highest_score)
        {
            second_highest = *highest_score;
            *highest_score = matchedp;
            highest_itr    = gmcitr;
        }else
        if(matchedp > second_highest)
        {
            second_highest = matchedp;
        }
        gmcitr ++;
    }
    //
    cout << "       : "    << this_raw_mkr        << "\n         with candidate insertion site between: "
         << (*highest_itr) << " and "             << *(highest_itr+1) << " with highest score " 
         << *highest_score << " (second-highest=" << second_highest   << ")." 
         << endl;
    //
    highest_itr ++;
    assert(highest_itr != this_gm_contigs.end());
    //
    return  *highest_itr; // will insert before ins_itr                            
}                                    
//
string find_best_map_insertion_site(bool                rely_on_paired,
                                    string              this_raw_mkr,
                                    string              this_raw_pat, 
                                    vector<string>      this_gm_contigs,
                                    map<string, string> contigMarkerPM,
                                    long*               highest_score)
{
    // position in a genetic map before which this raw marker will be inserted
    vector<string>::iterator highest_itr;
    // 1st checking: pair-end marker
    vector<string> targetinfo = split_string(this_raw_mkr, '\t');
    string paired_raw_mkr("");
    paired_raw_mkr += targetinfo[0];
    paired_raw_mkr += "\t";    
    if(this_raw_mkr.find("right")!=std::string::npos)
    {
        paired_raw_mkr += "left";
    }else
    if(this_raw_mkr.find("left")!=std::string::npos)    
    {
        paired_raw_mkr += "right";
    }else
    {
        // should never happen
        cout << "   Error: unexpected contig id " << this_raw_mkr << endl;
        return "unexpected.contig.id";
    }
    highest_itr = std::find(this_gm_contigs.begin(),
                            this_gm_contigs.end(),
                            paired_raw_mkr);
    /* if the paired-contig marker already in genetic map,
           insert this raw to its left or right
    */
    if(rely_on_paired && highest_itr != this_gm_contigs.end())
    {
        if(highest_itr == this_gm_contigs.begin())
        {
            /* this raw could not match better to the right of paired highest_itr
               otherwise, it could have been mapped here.
               So we put this raw on the left of the highest_itr
            */
            cout << "       : "    << this_raw_mkr    << "\n         with candidate insertion site before: "
                 << (*highest_itr) << " according to paired-marker in map. "
                 << endl;   
            *highest_score = this_raw_pat.size()*2; // 895 pollens * 2                     
            return  *highest_itr;                   // will insert before highest_itr                                      
        }else
        if(highest_itr == this_gm_contigs.end()-1)        
        {
            /* this raw could not match better to the left of paired highest_itr
               otherwise, it could have been mapped here.
               So we put this raw on the right of the highest_itr
            */
            cout << "       : "    << this_raw_mkr    << "\n         with candidate insertion site after--: "
                 << (*highest_itr) << " according to paired-marker in map. "
                 << endl;
            *highest_score = this_raw_pat.size()*2; // 895 pollens * 2
            highest_itr ++;                         // now it becomes this_gm_contigs.end()           
            return  *highest_itr;                   // will insert before highest_itr
        }else
        {
            string left_ctgmkr    = *(highest_itr - 1);
            string right_ctgmkr   = *(highest_itr + 1);  
            assert(contigMarkerPM.find(paired_raw_mkr) != contigMarkerPM.end());                         
            assert(contigMarkerPM.find(left_ctgmkr)    != contigMarkerPM.end());   
            assert(contigMarkerPM.find(right_ctgmkr)   != contigMarkerPM.end());      
            string highest_ctgpat = contigMarkerPM[paired_raw_mkr];         
            string left_ctgpat    = contigMarkerPM[left_ctgmkr];
            string right_ctgpat   = contigMarkerPM[right_ctgmkr]; 
            
            long matchedpleft  = 0;
            matchedpleft  += find_PM_pattern_match(this_raw_pat, left_ctgpat, true);            
            matchedpleft  += find_PM_pattern_match(this_raw_pat, highest_ctgpat, true);
            long matchedpright = 0;
            matchedpright += find_PM_pattern_match(this_raw_pat, highest_ctgpat, true);
            matchedpright += find_PM_pattern_match(this_raw_pat, right_ctgpat, true);    
            
            if(matchedpleft > matchedpright)
            {
                cout << "       : "      << this_raw_mkr << "\n         with candidate insertion site between: "
                     << *(highest_itr-1) << " and "      << *(highest_itr)   << " according to paired-marker in map. "
                     << endl; 
                *highest_score = this_raw_pat.size()*2; // 895 pollens * 2                     
                return  *highest_itr;                   // will insert before ins_itr                                      
            }else
            if(matchedpleft < matchedpright)
            {
                cout << "       : "      << this_raw_mkr << "\n         with candidate insertion site between: "
                     << *(highest_itr)   << " and "      << *(highest_itr+1) << " according to paired-marker in map. "
                     << endl; 
                *highest_score = this_raw_pat.size()*2; // 895 pollens * 2      
                highest_itr ++;               
                return  *highest_itr;                   // will insert before ins_itr                                      
            }else
            {
                cout << "       : a tie found for "        << this_raw_mkr 
                     << " on both side of paired marker "  << *(highest_itr) 
                     << ", inserted according to original contig left-right order. "
                     << endl; 
                if(paired_raw_mkr.find("left") != std::string::npos)
                {                    
                    highest_itr ++;
                }
                cout << "       : "        << this_raw_mkr << "\n         with candidate insertion site between: "
                     << *(highest_itr-1)   << " and "      << *(highest_itr) << " according to paired-marker in map. "
                     << endl;                                     
                *highest_score = this_raw_pat.size()*2;    // 895 pollens * 2                     
                return  *highest_itr;                      // will insert before ins_itr                    
            }        
        }
    }
    // 2nd checking: other markers
    vector<string>::iterator gmcitr;
    vector<string>::iterator gmcitr_end;
    gmcitr                       =    this_gm_contigs.begin();
    gmcitr_end                   = -- this_gm_contigs.end(); // last second 
    *highest_score  = 0;
    long second_highest = 0;
    highest_itr                  = gmcitr; // will insert between gmcitr and gmcitr+1
    while(gmcitr != gmcitr_end)
    {
        string left_ctgmkr  = *gmcitr;
        string right_ctgmkr = *(gmcitr+1);  
        //
        vector<string> leftinfo  = split_string(left_ctgmkr, '\t');
        vector<string> rightinfo = split_string(right_ctgmkr, '\t');
        if(leftinfo[0].compare(rightinfo[0])==0)
        {
            // do not insert between the same contig; namely 
            //    not between a\tright and b\tleft when a==b
            gmcitr ++;
            continue;
        }
        //
        assert(contigMarkerPM.find(left_ctgmkr)  != contigMarkerPM.end());   
        assert(contigMarkerPM.find(right_ctgmkr) != contigMarkerPM.end());               
        string left_ctgpat  = contigMarkerPM[left_ctgmkr];
        string right_ctgpat = contigMarkerPM[right_ctgmkr];
        //
        long matchedp = 0;
        matchedp += find_PM_pattern_match(this_raw_pat, left_ctgpat, true);
        matchedp += find_PM_pattern_match(this_raw_pat, right_ctgpat, true);
        if(matchedp > *highest_score)
        {
            second_highest = *highest_score;
            *highest_score = matchedp;
            highest_itr    = gmcitr;
        }else
        if(matchedp > second_highest)
        {
            second_highest = matchedp;
        }
        gmcitr ++;
    }
    //
    cout << "       : "    << this_raw_mkr        << "\n         with candidate insertion site between: "
         << (*highest_itr) << " and "             << *(highest_itr+1) << " with highest score " 
         << *highest_score << " (second-highest=" << second_highest   << ")." 
         << endl;
    //
    highest_itr ++;
    assert(highest_itr != this_gm_contigs.end());
    //
    return  *highest_itr; // will insert before ins_itr
}
long find_PM_pattern_match(string pm1, string pm2, bool scale)
{
    // pm1 is the target marker to be inserted in genetic map
    std::transform(pm1.begin(), pm1.end(), pm1.begin(), ::toupper);    
    std::transform(pm2.begin(), pm2.end(), pm2.begin(), ::toupper);    
    // only consider effective exact P/M positions
    long matchedp   = 0;
    unsigned long effectivep = 0;
    assert(pm1.size() == pm2.size());
    for(unsigned long ii = 0; ii < pm1.size(); ii ++)
    {
        string letter1 = pm1.substr(ii, 1);
        string letter2 = pm2.substr(ii, 1);
        //
        if(letter1.compare("U")!=0 && letter2.compare("U") != 0)
        {
            effectivep ++;
        } 
        //
        if(letter1.compare("P")==0 && letter2.compare("P") == 0)
        {
            matchedp += 2;
        }else
        if(letter1.compare("M")==0 && letter2.compare("M") == 0)
        {
            matchedp += 2;
        }else 
        if((letter1.compare("M")==0 && letter2.compare("P") == 0) ||
           (letter1.compare("P")==0 && letter2.compare("M") == 0) )
        {
            matchedp -= 1;
        }
        else;
    }
    if(effectivep==0) effectivep = 1; // avoid 0
    if(scale)
    {
        return (long) round( matchedp * 1.0 / effectivep * pm1.size() );
    }else
    {
        size_t n = std::count(pm1.begin(), pm1.end(), 'P');
        if(n == 0)
        {
            n = std::count(pm1.begin(), pm1.end(), 'M');
        }
        // in case the pattern in genetic map is with limited pollen information
        size_t m = std::count(pm2.begin(), pm2.end(), 'U');
        m = pm2.size() - m; // P or M
        if(n > m && m >= pm1.size()/4) // why 1/3? at least ~200 effective markers...
        {
            n = m;
        }
        if(n == 0) 
        {
            n = 1;
        }        
        return (long) round( matchedp * 1.0 / n  * pm1.size() );
    }
}
//
bool get_contig_PMsequence(string mkrfile, map<string, bool> ctgJoinMapPhasing_all, map<string, string>* contigMarkerPM)
{
    // returns contigMarkerPM with flipped PM patterns according to JoinMap phase info.
    ifstream ifp;
    ifp.open(mkrfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open PM pattern file " << mkrfile << endl;
        return false;
    }
    int phased_corrected = 0;    
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]  =='#') continue;
        if(line.find("left")==std::string::npos &&
           line.find("right")==std::string::npos)
        {
            continue;
        }
        //
        vector<string> lineinfo = split_string(line, '\t');
        string pmpat  = lineinfo[0];
        string contig = lineinfo[1];
        if(line.find("left") != std::string::npos)
        {
            contig += "\tleft";
            map<string, bool>::iterator phaseitr = ctgJoinMapPhasing_all.find(contig);            
            if(phaseitr != ctgJoinMapPhasing_all.end())
            {
                if((*phaseitr).second == true)
                {
                    pmpat = flip_PMpat(pmpat);
                    phased_corrected ++;
                }
            }
            (*contigMarkerPM).insert(std::pair<string, string>(contig, pmpat));
        }else
        if(line.find("right") != std::string::npos)
        {
            contig += "\tright";
            map<string, bool>::iterator phaseitr = ctgJoinMapPhasing_all.find(contig);            
            if(phaseitr != ctgJoinMapPhasing_all.end())
            {
                if((*phaseitr).second == true)
                {
                    pmpat = flip_PMpat(pmpat);
                    phased_corrected ++;
                }
            }            
            (*contigMarkerPM).insert(std::pair<string, string>(contig, pmpat));            
        }else ;
    } 
    cout << "   Info: " << (*contigMarkerPM).size() << " contig marker PM pattens collected; "   << endl; 
    cout << "         " << phased_corrected         << " flipped according to JoinMap phasing. " << endl;
    ifp.close();
    return true;
}            
//
string flip_PMpat(string PMpat)
{
    string flipped = "";
    for(int i=0; i < PMpat.size(); i++)
    {
        if(PMpat.substr(i, 1).compare("P")==0)
        {
            flipped += "M";
        }else
        if(PMpat.substr(i, 1).compare("M")==0)
        {
            flipped += "P";
        }else
        {
            flipped += "u"; // indicating flipped; otherwise "U"
        }
    }
    return flipped;
}  
//
bool get_del_like_PMsequence(string                           del_like_file,
                 map<string, map<unsigned long, DELMARKER> >* del_like_marker_PM)
{
    // map<string, map<unsigned long, DELMARKER> >: <contig, <del-sta, {ctg, sta, end, haphom, pmpat} > >
    ifstream ifp;
    ifp.open(del_like_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open del-like PM pattern file: " << del_like_file << endl;
        return false;        
    }
    cout << "   Info: PM will be re-determined with RPKM value - so it can differ from raw PM pattern. "  << endl;
    //
    int update_del_marker_definition = 0;    
    string example_update = "";
    unsigned long updated_size       = 0;
    //
    unsigned long del_like_num = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        
        if(line.find("rpkm") != std::string::npos)
        {
            // PUPPUUPUUUPUUP tig00000026_pilon   1   13743   hom_high    589/895 pmpat
            vector<string> lineinfo = split_string(line, '\t');
            if(lineinfo.size() < 7)
            {
                cout << "   Warning: skipped unexpected line with insufficient info at " << line << endl;
                continue;
            }
            del_like_num ++;
            DELMARKER this_node;            
            this_node.ctg    = lineinfo[1];
            this_node.sta    = strtoul(lineinfo[2].c_str(), NULL, 0);
            this_node.end    = strtoul(lineinfo[3].c_str(), NULL, 0);
            this_node.haphom = lineinfo[4];
            this_node.pmpat  = "";  // will be updated with rpkm value
            // update hap_low with cell support - new on 2020-02-11
            vector<string> cellcount = split_string(lineinfo[5], '/');
            double cellsupport       = atof(cellcount[0].c_str());
            double celltotal         = atof(cellcount[1].c_str());
            double cellsupportratio  = cellsupport / celltotal; // cutoff to switch hap_low to hom_low: 450/895=0.503
            /* with 895 cells, for hap_high we observed avg cell support of 298 => 450 ~ 298 + 298/2 = 447 */
            if(lineinfo[4].compare("hap_low")==0 && cellsupportratio>hap_low_to_hom_low_reset)
            {
                // reset del-marker state
                this_node.haphom = "hom_low";
                update_del_marker_definition ++;
                updated_size += (this_node.end - this_node.sta + 1);
                if(example_update.size()==0)
                {
                    example_update = line;
                }
            }
            //
            for(int ii = 0; ii < lineinfo[0].size(); ii ++)
            {
                // if rpkm>1
                if(lineinfo[0].substr(ii, 1).compare("0")!=0 && 
                   lineinfo[0].substr(ii, 1).compare("1")!=0
                )
                {
                    this_node.pmpat += "P";
                }else
                {
                    this_node.pmpat += "U";
                }
            }      
            if( (*del_like_marker_PM).find(this_node.ctg) == (*del_like_marker_PM).end() )
            {
                map<unsigned long, DELMARKER> this_del_marker;
                this_del_marker.insert(std::pair<unsigned long, DELMARKER>(this_node.sta, this_node));
                (*del_like_marker_PM).insert(std::pair<string, map<unsigned long, DELMARKER> >
                                               (this_node.ctg, this_del_marker));
            }else
            {
                (*del_like_marker_PM)[this_node.ctg].insert(std::pair<unsigned long, DELMARKER>
                                                                     (this_node.sta, this_node));
            }
        }
    }
    ifp.close();
    //
    cout << "   Info: " << del_like_num  << " del-like markers collected from " 
         << (*del_like_marker_PM).size() << " contigs."
         << endl;
    cout << "   Warning: " << update_del_marker_definition << " hap_low markers have been reset as hom_low "
         << "because they have more than 0.503 cell support, which should be less for confident del regions; "
         << endl
         << "            this is " << updated_size << " bp in size. " 
         << endl
         << "            for example: " << example_update << endl;
    // 
    return true;
}                    
//
bool get_contig_genetic_map(string               gmfile,
                            map<string, GNMAP>*  gmContig,
                            vector<string>*      gmContigOrder,
                            map<string, string>* unmappedContig)
{
    /* genetic map of a group of contig markers:
       gmfile        -- genetic map file with ordering of contigs
       gmContig      -- <contig-marker, <marker-map-details> >
       gmContigOrder -- initial order of contig markers in this group
       
       QUESTION: how to get order of paired-end markers? Minus/plus strand of a contig?
          Currently, unmapped ones always inserted as original left-right order,
                     left is always inserted before right.
          -- 20200102: check its similarity to surrounding pm patterns (of contigs in original genetic map by JoinMap).
    */
    ifstream ifp;
    ifp.open(gmfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gmfile << endl;
        return false;
    }
    vector<string> fileinfo = split_string(gmfile, '/');
    string groupLabel       = fileinfo[fileinfo.size()-1];
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        
        vector<string> lineinfo = split_string(line, ' ');
        if(lineinfo.size()<2)
        {
            cout << "   Warning: insufficient line info: " << line << endl;
            continue;
        }
        //
        string thiscontig = retrieve_contigid(lineinfo[0]);
        string leri       = "";
        if(lineinfo[0].find("_ri")!=std::string::npos)
        {
            leri          = "right";
        }else
        if(lineinfo[0].find("_le")!=std::string::npos)        
        {
            leri          = "left";
        }else
        {
            cout << "    Error: unexpected contig format: "            << line << endl;
            cout << "         : --expected contig format: tigxxxxxxxx_le/ri. " << endl;
            return false;
        }
        if(line.find("unmapped")!=std::string::npos)
        {
            (*unmappedContig).insert(std::pair<string, string>(thiscontig+"\t"+leri, "unmapped"));
            //cout << "   Check: " << thiscontig+"\t"+leri << " inserted in unmapped list. " << endl;
        }else
        {
            GNMAP tmpmap;
            tmpmap.distance = atof(lineinfo[1].c_str()); // original labeled genetic distance value 
            tmpmap.group    = groupLabel;                // linkage group
            tmpmap.status   = "JoinMapped";              // this is in map according to JoinMap
            // keep order of contigs in map
            (*gmContigOrder).push_back(thiscontig+"\t"+leri);
            (*gmContig).insert(std::pair<string, GNMAP>(thiscontig+"\t"+leri, tmpmap));
            // cout << "   Check: "    << thiscontig+"\t"+leri 
            //      << " inserted in   mapped list: gmval=" 
            //      << tmpmap.distance << ",  lineinfo=" << lineinfo.size() << endl;            
        }
    }
    cout << "   Info: in " << groupLabel << ", found " << (*unmappedContig).size() 
         << " unmapped contig markers;"  << endl;
    cout << "   Info: in " << groupLabel << ", found " << (*gmContigOrder).size()  
         << "   mapped contig markers."  << endl;
    if((*unmappedContig).size() > 0)
    {
        cout << "   Info: updating genetic maps based on paired-end marker info of a contig. " << endl;    
        map<string, string>::iterator unmapitr;
        map<string, string>::iterator unmapitr_end;
        unmapitr     = (*unmappedContig).begin();
        unmapitr_end = (*unmappedContig).end();
        while(unmapitr != unmapitr_end)
        {
            string unmappedcontig  = (*unmapitr).first;                  // unmapped contig \t left/right
            vector<string> ctginfo = split_string(unmappedcontig, '\t'); // unmapped contig \t left/right
            string pairendctg = ctginfo[0] + "\t";
            if(ctginfo[1].compare("right")==0)
            {
                pairendctg += "left";
            }else
            {
                pairendctg += "right";
            }
            vector<string>::iterator mapitr = std::find((*gmContigOrder).begin(), (*gmContigOrder).end(), pairendctg);
            if(mapitr != (*gmContigOrder).end())
            {
                //cout << "   Check: pair-end marker of " << pairendctg << " found in genetic map for unmapped: " 
                //     << unmappedcontig                  << endl;
                // update mapped list
                if(ctginfo[1].compare("right")==0)
                {
                    mapitr ++;                
                    (*gmContigOrder).insert(mapitr, unmappedcontig);
                    //cout << "   Check: " << unmappedcontig << " inserted after mapped " << pairendctg << endl;
                }else
                {
                    (*gmContigOrder).insert(mapitr, unmappedcontig);
                    //cout << "   Check: " << unmappedcontig << " inserted before mapped " << pairendctg << endl;                
                }
                // update mapped list with details of original unmapped contig
                GNMAP tmpmap;
                tmpmap.distance = (*gmContig)[pairendctg].distance; // original labeled genetic distance value 
                tmpmap.group    = (*gmContig)[pairendctg].group;    // linkage group
                tmpmap.status   = "pairLocated";                    // this is in map according to its paired-end marker
                // insert this unmapped contigs in map
                (*gmContig).insert(std::pair<string, GNMAP>(unmappedcontig, tmpmap));         
                //
                (*unmappedContig).erase(unmapitr ++);               
            }
            else
            {
                //cout << "   Check: " << unmappedcontig << " with pair " << pairendctg << " not in map." << endl;
                unmapitr ++;
            }
        }
        cout << "   Info: in "               << groupLabel                                                   << ", " 
             << (*unmappedContig).size()     << " unmapped contig markers after pair-end marker correction;" << endl;
        cout << "   Info: in "               << groupLabel                                                   << ", " 
             << (*gmContigOrder).size()      << " mapped   contig markers after pair-end marker correction." << endl;
        if((*unmappedContig).size() > 0)
        cout << "   Info: remaining unmapped means both end-markers are not in original genetic map. "       << endl
             << "         they would be inserted into the genetic map later based on PM similarity. "        << endl;
    }
    //
    ifp.close();             
    //
    return true;
}
//
string retrieve_contigid(string jmcontigid)          
{
    string thiscontig = "";
    //
    size_t postig = jmcontigid.find("tig");
    size_t posund = jmcontigid.find("_", postig);
    string tigid  = jmcontigid.substr(postig, posund-1-postig+1);
    thiscontig = tigid + "_pilon";
    // caution: special cases
    if(thiscontig.find("tigsA003826")!=std::string::npos)
    {
        thiscontig = "tig00003826_pilonsA";
        cout << "         special contig: " << thiscontig << endl;
    }else
    if(thiscontig.find("tig00003826")!=std::string::npos)
    {
        thiscontig = "tig00003826_pilonsB";
        cout << "         special contig: " << thiscontig << endl;
    }else
    if(thiscontig.find("tigsA003688")!=std::string::npos)
    {
        thiscontig = "tig00003688_pilonsA";
        cout << "         special contig: " << thiscontig << endl;                
    }else
    if(thiscontig.find("tig00003688")!=std::string::npos)
    {
        thiscontig = "tig00003688_pilonsB";
        cout << "         special contig: " << thiscontig << endl;
    }else ;
    //
    return thiscontig;
}
//
bool match_phasegroup_mapgroup(map<string, vector<string> >    contig_genetic_map, 
                               map<string, map<string, bool> > contig_phasing,
                               map<string, string>*            matched_phase_map,
                               map<string, string>*            matched_map_phase)
{
    /* match the group of phased contigs and the group of contigs in genetic map 
           as the files may in different orders from the input, e.g.,       
           phase file 1.txt may not be related to map_group1.txt, but to map_group6.txt
       Note: elements of contig_genetic_map should be a subset of contig_phasing       
    */    
    map<string, map<string, bool> >::iterator phasefitr;
    map<string, map<string, bool> >::iterator phasefitr_end;
    phasefitr     = contig_phasing.begin();
    phasefitr_end = contig_phasing.end();
    while(phasefitr != phasefitr_end)
    {
        string this_pg_file               = (*phasefitr).first;
        map<string, bool> this_pg_contigs = (*phasefitr).second;
        //
        map<string, string> matchfile;  // <phased_group_file, genetic_map_file>
        map<string, int>    matchscore; // <phased_group_file, matched_ctg_numb>                    
        //
        map<string, vector<string> >::iterator mapfitr;
        map<string, vector<string> >::iterator mapfitr_end;
        mapfitr     = contig_genetic_map.begin();
        mapfitr_end = contig_genetic_map.end();
        while(mapfitr != mapfitr_end)
        {
            string this_gmfile             = (*mapfitr).first; // only filename without path info
            vector<string> this_gm_contigs = (*mapfitr).second;
            //
            int this_score = 0; // how many contigs match         
            //
            vector<string>::iterator mapctgitr;
            vector<string>::iterator mapctgitr_end;            
            mapctgitr     = this_gm_contigs.begin();
            mapctgitr_end = this_gm_contigs.end();            
            while(mapctgitr != mapctgitr_end)
            {
                map<string, bool>::iterator phasectgitr = this_pg_contigs.find(*mapctgitr);
                if(phasectgitr != this_pg_contigs.end())
                {
                    this_score ++;
                }
                mapctgitr ++;
            } 
            if(this_score>0)
            cout << "   Check: matched "   
                 << this_score             << " between: phased group " 
                 << this_pg_file           << " (with "   
                 << this_pg_contigs.size() << " markers) and genetically ordered/mapped group " 
                 << this_gmfile            << " (with "   
                 << this_gm_contigs.size() << " markers)." << endl;  
            //
            if(matchfile.size()>0)
            {
                if(matchscore[this_pg_file]<this_score)
                {
                    matchfile.clear();
                    matchscore.clear();
                    matchfile.insert(std::pair<string, string>(this_pg_file, this_gmfile));
                    matchscore.insert(std::pair<string, int>(this_pg_file, this_score));                    
                }
            }else
            { 
                matchfile.insert(std::pair<string, string>(this_pg_file, this_gmfile));
                matchscore.insert(std::pair<string, int>(this_pg_file, this_score));
            }
            //                 
            mapfitr ++;
        }
        //
        (*matched_phase_map).insert(matchfile.begin(), matchfile.end());
        //
        phasefitr ++;        
    }    
    //
    map<string, string>::iterator mitr;
    map<string, string>::iterator mitr_end;
    mitr     = (*matched_phase_map).begin();
    mitr_end = (*matched_phase_map).end();
    while(mitr != mitr_end)
    {
        (*matched_map_phase).insert(std::pair<string, string>( (*mitr).second, (*mitr).first ) );        
        mitr ++;
    }
    //     
    return true;
}                 
//
bool get_contig_phase(string               jmfile, 
                      map<string, bool>*   ctgJoinMapPhasing,
                      map<string, bool>*   ctgJoinMapPhasing_all,
                      map<string, string>* ctgGroup,
                      int*                 nflip)
{
    // get ids with the respective phasing status of effective contigs in pollens
    ifstream jmifp;
    jmifp.open(jmfile.c_str(), ios::in);
    if(!jmifp.good())
    {
        cout << "   Error: cannot open file " << jmfile << endl;
        return false;
    }
    //
    vector<string> this_file_info = split_string(jmfile, '/');
    string jmflagstring = this_file_info[this_file_info.size()-1];
    /*
        414	tig00004185_ri	{0}	...
        413	tig00004185_le	{0}	...
        165	tig00003747_le	{1}	...
    */
    (*nflip) = 0;    
    while(jmifp.good())
    {
        string line("");
        getline(jmifp, line);
        if(line.size() == 0 ) continue;
        if(line[0]=='#')      continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<3) 
        {
            cout << "   warning: unexpected line: " << endl;
            continue;
        }
        //
        if(line.find("{0}")!=std::string::npos || line.find("{1}")!=std::string::npos)
        {
            size_t postig = line.find("tig");
            size_t posund = line.find("_", postig);
            
            string tigid  = line.substr(postig, posund-1-postig+1);
            
            string contigid = retrieve_contigid(tigid); // => tig00004185_pilon            
            
            string thiskey = "";
            bool   flip    = true;                        
            if(line.find("le")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = contigid + "\t" + "left";
                flip    = true;   
                (*nflip) ++;                 
            }else
            if(line.find("le")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = contigid + "\t" + "left";
                flip    = false;                    
            }else
            if(line.find("ri")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = contigid + "\t" + "right";
                flip    = true;   
                (*nflip) ++;                 
            }else
            if(line.find("ri")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = contigid + "\t" + "right";
                flip    = false;                    
            }else ;
            //
            // cout << "   Check: phasing contig key collected: " << thiskey << endl;
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
    //cout << "   Info: PM pattern of " << (*nflip) << " contig markers need flipping. " << endl;
    jmifp.close();
    //    
    return true;
}
//
bool get_group_files(string listfiles, vector<string>* groupfiles, string gmapphasing)
{
    ifstream ifp;
    ifp.open(listfiles.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << listfiles << endl;
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
            cout << "   Warning: file link not valid, skipped: " << line << endl;
            
        }else
        {
            (*groupfiles).push_back(line);
            cout << "   Info: " << gmapphasing << " file collected: " << line << endl;
            tmpifp.close();
        }
    }
    ifp.close();
    if((*groupfiles).size()==0) return false;
    else return true;
}
//
bool get_options(int                                     argc,
                 char*                                   argv[],
                 string*                                 fmap,
                 string*                                 fphase,
                 string*                                 fmarker,
                 string*                                 fmarker_del_like,
                 string*                                 outprefix, 
                 bool*                                   verbose)
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
        if(optstr.compare("--map") == 0)
        {
            ic ++;
            *fmap = (string)argv[ic];
            ifstream fp;
            fp.open((*fmap).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: genetic map file (=genetic ordering of contigs, e.g, by JoinMap) provided: "         
                     << *fmap << endl;
            }
            else
            {
                cout << "   Error: cannot open linkage map file "      
                     << *fmap << endl;
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
                cout << "   Info: contig marker file (=PM patterns at contig-paired-ends by snp markers) provided: "         
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
        if(optstr.compare("--marker-del-like") == 0)
        {
            ic ++;
            *fmarker_del_like = (string)argv[ic];
            ifstream fp;
            fp.open((*fmarker_del_like).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: del-like marker file (=PM patterns of regions without snp markers) provided: "         
                     << *fmarker_del_like << endl;
            }
            else
            {
                cout << "   Error: cannot open del-like marker file "      
                     << *fmarker_del_like << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--hap2hom") == 0)
        {
            ic ++;
            hap_low_to_hom_low_reset = atof(argv[ic]);
            cout << "   Info: hap_low region would be reset as hom_low if its cell support ratio is larger than "
                 << hap_low_to_hom_low_reset << endl;
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
                cout << "   Info: contig phasing file (=grouping of contigs according to contig-end-markers) provided: "         
                     << *fphase << endl;
            }
            else
            {
                cout << "   Error: cannot open contig phasing file "      
                     << *fphase << endl;
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
    cout << endl;
    // check necessary files - TODO
    return true;
}
