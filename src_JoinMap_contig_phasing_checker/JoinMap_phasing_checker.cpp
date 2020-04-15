/* given 2 genotype strings from asPollinator, i.e.,

	s2_genotype_contig_seq_parentalmarkers.txt
	s2_genotype_contig_seq.txt

   for the same set of pollen samples, 
   
   and 
   
   corresponding JoinMap contigs phasing status from Jose,
   
   	2019.10.17.pb_pollen_markers_0.1_nan_filtered_Group4_Map2_unicode.txt
   	(caution on window/linux file conversion: line ending might be different resulting in errors!)
   
   check the differences in genotype determination.
   
   Written by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
   
*/
#include       <map>
#include    <string>
#include    <vector>
#include   <fstream>
#include   <sstream>
#include  <iostream>
#include <algorithm>
#include  <string.h>
#include  <stdlib.h>
#include  <assert.h>
#include "split_string.h"

struct CONTIG
{
    string left;
    string co;
    string right;
    bool   paired; // whether the info exists in both sets
    bool   leftswap;
    bool   rightswap;
};
bool read_ctg_genotype(string               gtfile,   
                       map<string, CONTIG>* ctgGT);
int  get_genotype_diff(string               pattern1, 
                       string               pattern2, 
                       string*              match,
                       int*                 mismatch,
                       int*                 bothExist);                     
string swap_pattern(   string               pattern); // if JoinMap phasing status is {1}
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        cout << "\nFunction: given 2 file of s2_genotype_contig_seq.txt from currot/orangered and rojopasion markers, and ";
        cout << "\n          JoinMap grouping table with phasing status of contigs, ";
        cout << "\n          check the differences in defining the contig genotypes in pollens. " << endl;
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "Usage:  joinmap_phasing_checker s2_genotype_contig_seq_parental.txt s2_genotype_contig_seq_rojopasion.txt Group4_JoinMap_unicode.txt outflag\n" << endl;
        return 1;
    }
    // 
    cout << "   ----  check consistency between \n   " << argv[1] << "\n   "
                                                       << argv[2] << "\n   "
                                                       << argv[3] << endl; 
    // get ids with the respective phasing status of effective contigs in pollens
    string jmfile = (string)argv[3];
    ifstream jmifp;
    jmifp.open(jmfile.c_str(), ios::in);
    if(!jmifp.good())
    {
        cout << "   Error: cannot open file " << jmfile << endl;
        return 1;
    }
    /*
    (Individual:)			2	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	23	24	25	26	27	28	29	30	32	33	34	35	36	38	39	40	41	43	44	45	46	47	48	49	51	52	53	55	56	57	59	61	62	63	64	65	66	67	68	69	70	71	72	73	74	76	77	79	80	81	82	84	86	87	88	90	91	93	95	97	98	99	100	101	102	103	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	124	126	127	128	129	130	131	132	133	134	135	136	137	138	140	142	143	144	145	146	147	148	149	152	153	154	155	156	157	159	160	161	162	164	165	166	167	168	169	170	172	173	174	175	176	177	178	179	180	181	184	185	186	187	188	189	190	192	193	194	195	196	197	198	199	200	202	204	205	206	207	208	209	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	230	231	233	234	235	237	238	239	240	241	242	243	244	245	246	247	248	249	251	252	253	255	256	257	258	259	260	261	262	263	264	265	266	267	268	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	297	299	302	303	304	305	306	307	308	309	311	312	313	314	315	316	317	319	320	321	322	323	324	325	326	327	328	329	330	331	332	334	335	336	337	338	339	340	342	344	345	346	347	348	349	350	353	354	356	357	359	360	361	362	364	365	366	367	368	369	370	371	373	375	377	379	380	381	382	383	384	387	389	391	392	393	395	396	397	398	399	402	403	404	405	406	407	408	410	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	428	429	430	432	433	434	435	436	438	439	440	441	442	443	444	445
    414	tig00004185_ri	{0}	0	a	b	b	a	b	a	b	a	a	b	a	b	a	b	b	a	a	a	b	a	a	b	b	b	a	a	a	a	a	b	b	b	b	b	b	b	a	b	a	a	a	a	a	b	a	b	a	a	b	a	a	a	a	a	a	a	b	a	a	b	b	a	a	b	a	b	b	a	a	a	b	a	a	b	b	b	b	a	b	b	a	a	a	a	a	b	a	b	b	a	a	a	b	a	a	b	b	a	a	b	a	b	b	a	a	b	a	b	a	b	b	b	b	a	a	a	b	b	a	b	b	a	b	a	b	a	a	b	a	a	a	b	a	a	a	a	b	b	a	a	a	a	a	b	b	b	b	b	a	b	b	a	a	a	b	a	a	a	a	b	a	b	b	a	b	a	a	b	a	a	b	a	a	a	a	a	b	b	b	b	b	a	b	a	b	b	b	b	b	b	a	b	a	a	b	b	a	a	b	b	a	b	a	b	a	b	a	a	a	a	a	a	a	a	a	b	b	b	b	b	b	a	a	a	a	a	-	a	a	a	b	a	a	a	b	b	b	a	b	a	a	b	a	a	a	b	b	b	a	a	a	a	a	b	a	b	b	b	b	a	a	b	a	a	b	a	a	b	b	b	a	a	a	-	b	a	a	a	a	b	a	a	b	b	b	b	b	a	b	a	a	a	a	b	b	b	b	a	b	a	b	b	a	a	b	a	a	b	a	b	a	a	a	b	a	a	a	b	b	a	a	b	a	b	b	a	b	a	a	b	a	a	b	b	b	a	b	a	a	b	b	a	b	b	a	a	b	a	b	a	b	b	a	a	b	b	b	b	b	a	b	a	b	a	b	b	b	a	b	a	b	b	b	a	b	b
    413	tig00004185_le	{0}	0	a	b	b	a	b	a	b	a	a	b	a	b	a	b	b	a	a	a	b	a	a	b	b	b	a	a	a	a	a	b	b	b	b	b	b	b	a	b	a	a	a	a	a	b	a	b	a	a	b	a	a	a	a	a	a	a	b	a	a	b	b	a	a	b	a	b	b	a	a	a	b	a	a	b	b	b	b	a	b	b	a	a	a	a	a	b	a	b	b	a	a	a	b	a	a	b	b	a	a	b	a	b	b	a	a	b	a	b	a	b	b	b	b	a	a	a	b	b	a	b	b	a	b	a	b	a	a	b	a	a	a	b	a	a	a	a	b	b	a	a	a	a	a	b	b	b	b	b	a	b	b	a	a	a	b	a	a	a	a	b	a	b	b	a	b	a	a	b	a	a	b	a	a	a	a	a	b	b	b	b	b	a	b	a	b	b	b	b	b	b	a	b	a	a	b	b	a	a	b	b	a	b	a	b	a	b	a	a	a	a	a	a	a	a	a	b	b	b	b	b	b	a	a	a	a	a	-	a	a	a	b	a	a	a	b	b	b	a	b	a	a	b	a	a	a	b	b	b	a	a	a	a	a	b	a	b	b	b	b	a	a	b	a	a	b	a	a	b	b	b	a	a	a	-	b	a	a	a	a	b	a	a	b	b	b	b	b	a	b	a	a	a	a	b	b	b	b	a	b	a	b	b	a	a	b	a	a	b	a	b	a	a	a	b	a	a	a	b	b	a	a	b	a	b	b	a	b	a	a	b	a	a	b	b	b	a	b	a	a	b	b	a	b	b	a	a	b	a	b	a	b	b	a	a	b	b	b	b	b	a	b	a	b	a	b	b	b	a	b	a	b	b	b	a	b	b
    165	tig00003747_le	{1}	0.013	b	a	a	b	a	b	a	b	b	a	b	a	b	a	a	b	b	b	a	b	b	a	a	a	b	b	b	b	b	a	a	a	a	a	a	a	b	a	b	b	b	b	b	a	b	a	b	b	a	b	b	b	b	b	b	b	a	b	b	a	a	b	b	a	b	a	a	b	b	b	a	b	b	a	a	a	a	b	a	a	b	b	b	b	b	a	b	a	a	b	b	b	a	b	b	a	a	b	b	a	b	a	a	b	b	a	b	a	b	a	a	a	a	b	b	b	a	a	b	a	a	b	a	b	a	b	b	a	b	b	b	a	b	b	b	b	a	a	b	b	b	b	b	a	a	a	a	a	b	a	a	b	b	b	a	b	b	b	b	a	b	a	a	b	a	b	b	a	b	b	a	b	b	b	b	b	a	a	a	a	a	b	a	b	a	a	a	a	a	a	b	a	b	b	a	a	b	b	a	a	b	a	b	a	b	a	b	b	b	b	b	b	b	b	b	a	a	a	a	a	a	b	b	b	b	b	a	b	b	b	a	b	b	b	a	a	a	b	a	b	b	a	b	b	b	a	a	a	b	b	b	b	b	a	b	a	a	a	a	b	b	a	b	b	a	b	b	a	a	a	b	b	b	a	a	b	b	b	b	a	b	b	a	a	a	a	a	b	a	b	b	b	b	a	a	a	a	b	a	b	a	a	b	b	a	b	b	a	b	a	b	b	b	a	b	b	b	a	a	b	b	a	b	a	a	b	a	b	b	a	b	b	a	a	a	b	a	b	b	a	a	b	a	a	b	b	a	b	a	b	a	a	b	b	a	a	a	a	a	b	a	b	a	b	a	a	a	b	a	b	a	a	a	b	a	a
    */ 
    map<int, int>     orderpid;           // map<column-interger-id, pollen-interger-id>
    map<string, bool> ctgJoinMapPhasing;
    vector<string>    ordercid;           // vector<ctgid#left/right>, order keeped
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
        if(lineinfo[0].compare("(Individual:)")==0)
        {
            cout << "   Info: initialize pollen ids: ";
            for(int pid=1; pid<lineinfo.size(); pid ++)
            {
                int realpid = atoi(lineinfo[pid].c_str());
                orderpid.insert(std::pair<int, int>(pid-1, realpid));
            }
            cout << endl;
        }
        if(line.find("{0}")!=std::string::npos || line.find("{1}")!=std::string::npos)
        {
            size_t postig = line.find("tig");
            size_t posund = line.find("_", postig);
            
            string tigid  = line.substr(postig, posund-1-postig+1);
            
            string thiskey = "";
            bool   flip    = true;                        
            if(line.find("le")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon#left";
                flip    = true;                    
            }else
            if(line.find("le")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon#left";
                flip    = false;                    
            }else
            if(line.find("ri")!=std::string::npos && line.find("{1}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon#right";
                flip    = true;                    
            }else
            if(line.find("ri")!=std::string::npos && line.find("{0}")!=std::string::npos)
            {
                thiskey = tigid + "_pilon#right";
                flip    = false;                    
            }else ;
            //
            assert(ctgJoinMapPhasing.find(thiskey) == ctgJoinMapPhasing.end());
            ctgJoinMapPhasing.insert(std::pair<string, bool>(thiskey, flip));     
            // cout << "   check: " << thiskey << " as " << flip << endl;            
            //
            ordercid.push_back(thiskey);
        }
    }
    jmifp.close();
    //
    string outflag = (string)argv[4];
    // check parent input
    string parent = (string)argv[1];   
    ifstream pifp;
    pifp.open(parent.c_str(), ios::in);
    if(!pifp.good())
    {
        cout << "   Error: cannot open file " << parent << endl;
        return 1;
    }
    pifp.close();
    // check child input
    string child = (string)argv[2];
    ifstream rifp;
    rifp.open(child.c_str(), ios::in);
    if(!rifp.good())
    {
        cout << "   Error: cannot open file " << child << endl;
        return 1;
    }
    rifp.close();       
    // prepare paired output file
    vector<string> grpfileinfo = split_string(jmfile, '/');
    string outfile = "res_checking_parent_child_contig_phasing_" + grpfileinfo[grpfileinfo.size()-1];
    ofstream ofp;
    ofp.open(outfile.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file to write result: " << outfile << endl;
        return 1;
    }
    // step 1. read contig genotype
    map<string, CONTIG> ctgGTparent;
    if(!read_ctg_genotype(parent, &ctgGTparent) )
    {
       cout << "   Error: failed in reading " << parent << endl;
       return 1;
    }
    map<string, CONTIG> ctgGTchild;
    if(!read_ctg_genotype(child, &ctgGTchild) )
    {
       cout << "   Error: failed in reading " << child << endl;
       return 1;
    }    
    // step 2. compare two sets
    map<string, string> selected;       // <ctg#left/right, info-output>
    map<string, CONTIG>::iterator citr; // child itr
    map<string, CONTIG>::iterator citr_end;
    citr     = ctgGTchild.begin();
    citr_end = ctgGTchild.end();
    while(citr != citr_end)
    {
        string ctgID     = (*citr).first;
        CONTIG childctg  = (*citr).second;
        map<string, CONTIG>::iterator pitr = ctgGTparent.find(ctgID); // parent itr
        //
        if(pitr != ctgGTparent.end() ) // paired
        {
            (*citr).second.paired = true;
            (*pitr).second.paired = true;
            //
            // 2.1 comparing left sequence
            string thiskey  = ctgID + "#left";
            if(ctgJoinMapPhasing.find(thiskey) != ctgJoinMapPhasing.end())
            {
                int childeff = 0;  // effective child contig-markers with clear P/M, but not U
                childeff    += std::count((*citr).second.left.begin(), (*citr).second.left.end(), 'P');
                childeff    += std::count((*citr).second.left.begin(), (*citr).second.left.end(), 'M');
                // left sequence (possibly with child pattern swapped)
                string leftpat  = (*citr).second.left;
                string thisflip = "nonflipped";
                if(ctgJoinMapPhasing[thiskey] == true)
                {
                    leftpat  = swap_pattern(leftpat);
                    thisflip = "flipped";
                }
                string match2 = "";
                int mismatch2 = 0;
                int bothExist = 0; 
                int score2    = get_genotype_diff(leftpat, (*pitr).second.left,  &match2, &mismatch2, &bothExist);
                int parenteff = 0; // effective parent contig-markers with clear P/M, but not U
                parenteff += std::count((*pitr).second.left.begin(), (*pitr).second.left.end(), 'P');
                parenteff += std::count((*pitr).second.left.begin(), (*pitr).second.left.end(), 'M');
                //
                std::stringstream ss;
                ss.str("");
                ss << leftpat             << "\t" << (*citr).first  << "\tleft-sequence" << "\tchild\t"  << childeff  << "\t" << thisflip                << endl;
                ss << match2              << "\t" << (*citr).first  << "\tleft-sequence" << "\tmisma\t"  << mismatch2 << "\t" << mismatch2*1.0/bothExist << endl;               
                ss << (*pitr).second.left << "\t" << (*pitr).first  << "\tleft-sequence" << "\tparent\t" << parenteff << "\t" << bothExist               << endl;                
                ss << endl;
                selected.insert(std::pair<string, string>(thiskey, ss.str()));
            }else
            {
                ; // not in JoinMap-phased group
            }
            // 2.3 comparing right sequence unchanged
            thiskey  = ctgID + "#right";
            if(ctgJoinMapPhasing.find(thiskey) != ctgJoinMapPhasing.end())            
            {            
                int childeff = 0;
                childeff += std::count((*citr).second.right.begin(), (*citr).second.right.end(), 'P');
                childeff += std::count((*citr).second.right.begin(), (*citr).second.right.end(), 'M');            
                // 2.4 comparing right sequence with child pattern swapped
                string rightpat = (*citr).second.right; 
                string thisflip = "nonflipped";                
                if(ctgJoinMapPhasing[thiskey] == true)
                {
                    rightpat = swap_pattern(rightpat);
                    thisflip = "flipped";                    
                }                  
                string match2 = "";
                int mismatch2 = 0;    
                int bothExist = 0;                 
                int score2    = get_genotype_diff(rightpat, (*pitr).second.right,  &match2, &mismatch2, &bothExist); //      
                int parenteff = 0;
                parenteff += std::count((*pitr).second.right.begin(), (*pitr).second.right.end(), 'P');
                parenteff += std::count((*pitr).second.right.begin(), (*pitr).second.right.end(), 'M');                
                //
                std::stringstream ss;
                ss.str("");                
                ss << rightpat             << "\t" << (*citr).first  << "\tright-sequence" << "\tchild\t"  << childeff  << "\t" << thisflip                << endl;
                ss << match2               << "\t" << (*citr).first  << "\tright-sequence" << "\tmisma\t"  << mismatch2 << "\t" << mismatch2*1.0/bothExist << endl;
                ss << (*pitr).second.right << "\t" << (*pitr).first  << "\tright-sequence" << "\tparent\t" << parenteff << "\t" << bothExist               << endl;                    
                //
                ss << endl;
                selected.insert(std::pair<string, string>(thiskey, ss.str()));                
            }else
            {
                ; // not in JoinMap-phased group
            }
        }else
        {
            // parental genotyping on the related contig marker is missing from inputs!
            (*citr).second.paired = false;
            //
            // 2.1 comparing left sequence
            string thiskey  = ctgID + "#left";
            cout << "   Info: this contig " << ctgID << " not found in parent genotype list. " << endl;                            
            if(ctgJoinMapPhasing.find(thiskey) != ctgJoinMapPhasing.end())
            {
                cout << "         but found in JoinMap phasing left-list. " << endl;         // in JoinMap-phased group                
                int childeff = 0;  // effective child contig-markers with clear P/M, but not U
                childeff    += std::count((*citr).second.left.begin(), (*citr).second.left.end(), 'P');
                childeff    += std::count((*citr).second.left.begin(), (*citr).second.left.end(), 'M');
                // left sequence (possibly with child pattern swapped)
                string leftpat  = (*citr).second.left;
                string thisflip = "nonflipped";
                if(ctgJoinMapPhasing[thiskey] == true)
                {
                    leftpat  = swap_pattern(leftpat);
                    thisflip = "flipped";
                }
                string match2 = "";
                int mismatch2 = 0;
                int bothExist = 0;                 
                int score2    = get_genotype_diff(leftpat, leftpat,  &match2, &mismatch2, &bothExist); //    
                //
                std::stringstream ss;
                ss.str("");
                ss << leftpat             << "\t" << (*citr).first  << "\tleft-sequence" << "\tchild\t"          << childeff  << "\t" << thisflip                << endl;
                ss << match2              << "\t" << (*citr).first  << "\tleft-sequence" << "\tmisma\t"          << mismatch2 << "\t" << mismatch2*1.0/bothExist << endl;                
                ss << leftpat             << "\t" << (*citr).first  << "\tleft-sequence" << "\tparent-missing\t" << childeff << "\t"  << bothExist               << endl;                
                ss << endl;
                selected.insert(std::pair<string, string>(thiskey, ss.str()));
            }else
            {
                cout << "         also not in found in JoinMap phasing left-list. " << endl; // not in JoinMap-phased group
            }
            // 2.3 comparing right sequence unchanged
            thiskey  = ctgID + "#right";
            if(ctgJoinMapPhasing.find(thiskey) != ctgJoinMapPhasing.end())            
            {            
                cout << "         but found in JoinMap phasing right-list. " << endl;        // in JoinMap-phased group                                
                int childeff = 0;
                childeff += std::count((*citr).second.right.begin(), (*citr).second.right.end(), 'P');
                childeff += std::count((*citr).second.right.begin(), (*citr).second.right.end(), 'M');            
                // 2.4 comparing right sequence with child pattern swapped
                string rightpat = (*citr).second.right;
                string thisflip = "nonflipped";        
                if(ctgJoinMapPhasing[thiskey] == true)
                {
                    rightpat = swap_pattern(rightpat);
                    thisflip = "flipped";    
                }
                string match2 = "";
                int mismatch2 = 0;
                int bothExist = 0;                 
                int score2    = get_genotype_diff(rightpat, rightpat,  &match2, &mismatch2, &bothExist); //
                //
                std::stringstream ss;
                ss.str("");
                ss << rightpat             << "\t" << (*citr).first  << "\tright-sequence" << "\tchild\t"          << childeff  << "\t" << thisflip                << endl;
                ss << match2               << "\t" << (*citr).first  << "\tright-sequence" << "\tmisma\t"          << mismatch2 << "\t" << mismatch2*1.0/bothExist << endl;
                ss << rightpat             << "\t" << (*citr).first  << "\tright-sequence" << "\tparent-missing\t" << childeff  << "\t" << bothExist               << endl;                    
                //
                ss << endl;
                selected.insert(std::pair<string, string>(thiskey, ss.str()));
            }else
            {
                cout << "         also not in found in JoinMap phasing right-list. " << endl; // not in JoinMap-phased group
            }               
        }
        citr ++;
    }
    // output selected contig info in the orignal input order
    vector<string>::iterator cinitr;
    vector<string>::iterator cinitr_end;
    cinitr     = ordercid.begin();
    cinitr_end = ordercid.end();
    while(cinitr != cinitr_end)
    {
        map<string, string>::iterator sitr = selected.find(*cinitr);
        if(sitr == selected.end())
        {
            cout << "   warning: " << (*sitr).first << " not found in selected. " << endl;
            cinitr ++;            
            continue;
        }
        ofp << (*sitr).second;
        cinitr ++;
    }    
    //
    return 0;
}
//
string swap_pattern(string pattern)
{
    string tmp = "";
    
    for(int ii = 0; ii < pattern.size(); ii ++)
    {
        string letter = pattern.substr(ii, 1);
        if(letter.compare("P")==0 || letter.compare("p")==0)
        {
            tmp += "M";
        }
        else
        if(letter.compare("M")==0 || letter.compare("m")==0)
        {
            tmp += "P";
        }else
        {
            tmp += "U";
        }        
    }
    return tmp;
}
//
int get_genotype_diff(string pattern1, string pattern2, string* match, int* mismatch, int* bothExist)
{
    int score = 0; 
    *mismatch = 0;
    *bothExist= 0;    
    (*match)  = "";   
    assert(pattern1.size() == pattern2.size());    
    for(int ii=0; ii < pattern1.size(); ii ++)
    {
        string letter1 = pattern1.substr(ii, 1);
        string letter2 = pattern2.substr(ii, 1);
        
        if(letter1.compare("U")!=0 && letter2.compare("U")==0) // in child not in parent
        {
            (*match) += "+";
        }else
        if(letter1.compare("U")==0 && letter2.compare("U")!=0) // not in child in parent
        {
            (*match) += "-";
        }else        
        if(letter1.compare(letter2) == 0)                      // same pattern
        {
            score ++;
            (*match) += "|";
        }else
        {
            (*match) += "*";           
        }
        //
        if( (letter1.compare("P")==0 && letter2.compare("M")==0) ||
            (letter1.compare("M")==0 && letter2.compare("P")==0) )
        {
            (*mismatch) ++;
        }    
        // markers where both parent and child give clear genotype P/M but not U
        if(letter1.compare("U")!=0 && letter2.compare("U")!=0)
        {
            (*bothExist) ++;
        }        
    }
    return score;
}
//
bool read_ctg_genotype(string gtfile, map<string, CONTIG>* ctgGT)
{
    ifstream ifp;
    ifp.open(gtfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gtfile << endl;
        return 1;
    }    
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]=='#')   continue;
        //
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<3)
        {
            cout << "   Warning: unexpected line " << line << endl;
            continue;
        }
        //
        string pattern  = lineinfo[0]; // PM pattern
        string contigid = lineinfo[1]; // contig id
        string lr_co    = lineinfo[2]; // left/right/co flag
        //
        map<string, CONTIG>::iterator ctgitr;
        ctgitr = (*ctgGT).find(contigid);
        //
        if(ctgitr == (*ctgGT).end())
        {
            CONTIG tmp;
            assert(lr_co.find("left-sequence")!=std::string::npos);
            tmp.left  = pattern;
            tmp.right = "";
            tmp.co    = "";
            tmp.paired= false;
            (*ctgGT).insert(std::pair<string, CONTIG>(contigid, tmp));
        }else
        {
            if(lr_co.find("right-sequence")!=std::string::npos)
            {
                (*ctgitr).second.right = pattern;
            }else
            if(lr_co.find("CO-info")!=std::string::npos)
            {
                (*ctgitr).second.co    = pattern;
            }else ;        
        }
    }
    cout << "   Info: " << (*ctgGT).size() << " contig genotype pattern collected. " << endl;
    ifp.close();
    return true;
}










