/*
    Given 
    
        two PM patterns, check how much of them match and mistach.
        
    2019-12-30 by Hequan Sun
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
//
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        // g++ match_reporter.cpp split_string.cpp -O3 -o match_reporter
        cout << "\n   Function: check how two PM patterns match and mistach. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";        
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "match_reporter PM.txt" << endl;
        cout << " PM.txt gives paired PM patterns to compare, every two lines make a pair. " << endl << endl;
        return 1;
    }
    double startT= clock();
   
    string ifile = argv[1];
    ifstream ifp;
    ifp.open(ifile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << ifile << endl;
        return 1;
    }
    while(ifp.good())
    {
        string line1("");
        getline(ifp, line1);
        if(line1.size()==0) continue;
        if(line1[0] == '#') continue;        
        string line2("");
        getline(ifp, line2);
        if(line2.size()==0 || line2[0] == '#')
        {
            cout << "   Error: incorrect file format. " << endl;
            cout << "   PM11 PM11 info" << endl;
            cout << "   PM12 PM12 info" << endl;
            cout << "   ..."            << endl;
            return 1;
        }
        vector<string> line1info = split_string(line1, '\t');
        vector<string> line2info = split_string(line2, '\t');        
          
        string PM1 = line1info[0];
        string PM2 = line2info[0];
        
        if(PM1.size() != PM2.size())
        {
            cout << "   Error: size of two patterns do not match. " << endl;
            return 1;
        }
        int mm = 0;
        int pp = 0;
        int pm = 0;
        int mp = 0;
        int ot = 0; // other
        string tmp = "";
        for(int ii = 0; ii < PM1.size(); ii ++)
        {
            if(PM1.substr(ii, 1).compare("M") == 0 && 
               PM2.substr(ii, 1).compare("M") == 0)
            {
                mm ++;
                tmp += "|";
            }            
            if(PM1.substr(ii, 1).compare("P") == 0 && 
               PM2.substr(ii, 1).compare("P") == 0)
            {
                pp ++;
                tmp += "|";
            }else
            if(PM1.substr(ii, 1).compare("P") == 0 && 
               PM2.substr(ii, 1).compare("M") == 0)
            {
                pm ++;
                tmp += "*";                
            }else
            if(PM1.substr(ii, 1).compare("M") == 0 && 
               PM2.substr(ii, 1).compare("P") == 0)
            {
                mp ++;
                tmp += "*";                                
            }else
            {
                ot ++;
                tmp += ".";   // "U" related             
            }                    
        }
        cout << PM1 << endl 
             << PM2 << endl;
        cout << "report: mm+pp+pm+mp+ot=" << mm+pp+pm+mp+ot << " vs pattern_len=" << PM1.size() << endl
             << "   mm="   << mm << endl
             << "   pp="   << pp << endl
             << "   pm="   << pm << endl
             << "   mp="   << mp << endl
             << "   ot="   << ot << endl << endl;                          
    }
    ifp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
