/* given    
        JoinMap phasing status of a group contigs from: 1.txt
   	(caution on window/linux file conversion: line ending might be different resulting in errors!)
   	
   and 
        a path to ./asPollinator-run/ before 
   update the PM values in /s6_PM_pollen_bed_ctgwise/s6_PM_region_pollens_at_contig_*_pilon.txt
       caution on 3688 and 3826 contigs which each was manually separated as two sub-contigs (ids-need-change).
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
#include     <iomanip>
#include  <sys/stat.h>
#include    <dirent.h>
#include "split_string.h"
bool swap_pattern(bool this_flip, string* PM, string* precluster);
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\nFunction: given JoinMap grouping table with phasing status of contigs, and the path to asPollinator-run";
        cout << "\n          correct the values of P/M in related /s6_PM_pollen_bed_ctgwise/s6_PM_region_pollens_at_contig_*_pilon.txt. " << endl;
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "Usage:  pollen_bed_PMflipper Group4_JoinMap_unicode.txt path/to/asPollinator-run/ grpid" << endl;
        cout << "        path is before at upper level of /s6_PM_pollen_bed_ctgwise/\n"                     << endl;
        return 1;
    }
    //
    cout << endl;
    cout << "   Info: correction of PM values with JoinMap contigs phasings started.. " << endl;    
    cout << "   ----  you always have to check consistency between \n   " << argv[1] << "\n      and\n   "
                                                                          << argv[2] << endl;                                                                         
    // get ids with the respective phasing status of effective contigs in pollens
    string jmfile = (string)argv[1];
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
    vector<string>    ordercid;           // vector<ctgid>, order keeped; is there different in phasing in "#left/right"? caution!!!! currently says no.
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
            //cout << "   Info: initialize pollen ids: ";
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
                cout << "   Info: special contig: " << thiskey << endl;
            }else
            if(thiskey.find("tig00003826")!=std::string::npos)
            {
                thiskey = "tig00003826_pilonsB";
                cout << "   Info: special contig: " << thiskey << endl;                
            }else            
            if(thiskey.find("tigsA003688")!=std::string::npos)
            {
                thiskey = "tig00003688_pilonsA";
                cout << "   Info: special contig: " << thiskey << endl;                
            }else
            if(thiskey.find("tig00003688")!=std::string::npos)
            {
                thiskey = "tig00003688_pilonsB";
                cout << "   Info: special contig: " << thiskey << endl;                
            }else ;
            //
            map<string, bool>::iterator tmpcitr = ctgJoinMapPhasing.find(thiskey);
            if(tmpcitr != ctgJoinMapPhasing.end())
            {
                if((*tmpcitr).second != flip)
                {
                    cout << "   Warning: left and right phasings differ on the same contig " << thiskey << endl;
                    cout << "            need to update my code!! "                          << endl;
                    // this should only happens on mis-assembled contigs?
                }
            }else
            {
                ctgJoinMapPhasing.insert(std::pair<string, bool>(thiskey, flip));    
                // cout << "   check: " << thiskey << " as " << flip << endl;
            }            
            //
            ordercid.push_back(thiskey);
        }
    }
    jmifp.close();
    //
    // step 2. correct PM settings of invidual contigs in the current group
    string workfolder = (string)argv[2];
    string groupid    = (string)argv[3];

    // subfolder for collecting bed information of P/M of pollens along a contig 
    string tmpfolder_s8 = workfolder+"/s6_PM_pollen_bed_ctgwise_JoinMap_flipped/";
    DIR* dir = opendir(tmpfolder_s8.c_str());
    if (dir)
    {
        /* Directory exists. */
        cout << "   Info: output folder exists, not changed. " << endl;
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder_s8.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder_s8 << endl;
            return false;
        }else
        {   
            cout << "   Info: folder created: " << tmpfolder_s8 << endl;
        }
    }
    else;      
    map<string, bool>::iterator phasedctgitr;
    map<string, bool>::iterator phasedctgitr_end;
    phasedctgitr     = ctgJoinMapPhasing.begin();
    phasedctgitr_end = ctgJoinMapPhasing.end();
    int updatefilenum = 0;
    while(phasedctgitr != phasedctgitr_end)
    {       
        string this_contig  = (*phasedctgitr).first; // format: tig00000013_pilon.txt
        bool   this_flip    = (*phasedctgitr).second;
        //
        string this_bedfile = workfolder + "/s6_PM_pollen_bed_ctgwise/s6_PM_region_pollens_at_contig_" + this_contig + ".txt";
        ifstream ifp;
        ifp.open(this_bedfile.c_str(), ios::in);
        if(!ifp.good())
        {
            cout << "   Error: in bed file not found: " << this_bedfile << endl;
            return 1;
        }
        //
        string upda_bedfile = workfolder + "/s6_PM_pollen_bed_ctgwise_JoinMap_flipped/s6_PM_region_pollens_at_contig_" + this_contig + "_grp"+groupid+".txt";
        ofstream ofp;
        ofp.open(upda_bedfile.c_str(), ios::out);
        if(!ofp.good())
        {
            cout << "   Error: out bed file not found " << upda_bedfile << endl;
            return 1;
        }
        //
        while(ifp.good())
        {
            string line("");
            getline(ifp, line);
            if(line.size()== 0 ) continue;
            //
            if(line[0]    =='#')
            {
                ofp << line << endl;
                continue;
            }else
            {
                // tig00000013_pilon	1	50728	1	56	M	0_m	4_x
                vector<string> lineinfo = split_string(line, '\t');
                string PM         = lineinfo[5];
                string precluster = lineinfo[6];
                if(!swap_pattern(this_flip, &PM, &precluster))
                {
                    cout << "   Error: at " << line << endl;
                    return 1;
                }
                for(int ii = 0; ii < lineinfo.size(); ii ++)
                {
                    if(ii == 5)
                    {
                        ofp << PM << "\t";
                    }else
                    if(ii == 6)
                    {
                        ofp << precluster << "\t";
                    }else
                    if(ii<lineinfo.size()-1)
                    {
                        ofp << lineinfo[ii] << "\t";
                    }else
                    {
                        ofp << lineinfo[ii] << endl;
                    }
                }
            }         
        }
        //
        updatefilenum ++;
        //
        ifp.close();
        ofp.close();
        //
        phasedctgitr ++;
    }
    //
    cout << "   Info: correction of PM values along " << updatefilenum
         << " contigs with JoinMap phasings successfully done. "   << endl << endl;
    //
    return 0;
}
//
bool swap_pattern(bool this_flip, string* PM, string* precluster)
{    
    if(this_flip==true)
    {
        if((*PM).compare("M")==0)
        {
            (*PM) = "P";
        }else
        if((*PM).compare("P")==0)   
        {
            (*PM) = "M";
        }else
        {
            cout << "   Error: unexpected PM-value: " << (*PM) << endl;
            return false;
        }
        //
        if((*precluster).compare("0_m")==0)
        {
            (*precluster) = "0_p";
        }else
        if((*precluster).compare("0_p")==0)  
        {
            (*precluster) = "0_m";
        }else
        {
            ;
            // "0_u exist!!!"
            // cout << "   Error: unexpected cluster-flag: " << (*precluster) << endl;
            // return false;            
        }
    }    
    return true;
}





