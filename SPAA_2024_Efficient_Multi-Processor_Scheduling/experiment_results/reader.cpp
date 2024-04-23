/*
Copyright 2024 Huawei Technologies Co., Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@author Pal Andras Papp
*/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;

const set<string> tiny = {"CG_N2_", "CG_N3_", "CG_N4_", "exp_N4_", "exp_N5_", "exp_N6_", "kNN_N4_", "kNN_N5_", "kNN_N6_", "spmv_N10_", "spmv_N6_", "spmv_N7_", "bicgstab", "NN_3_gyro", "means", "pregel_p"};
const set<string> small = {"CG_N5_", "CG_N6_", "CG_N7_", "CG_N8_", "CG_N9_K5_", "exp_N10_", "exp_N15_", "exp_N18_", "exp_N20_", "exp_N25_", "kNN_N10_", "kNN_N13_", "kNN_N15_", "kNN_N20_", "kNN_N25_", "spmv_N25_", "spmv_N35_", "spmv_N40_", "pregel_cc", "pagerank", "snni"};
const set<string> medium = {"CG_N12_", "CG_N15_", "CG_N17_", "CG_N21_", "CG_N9_K9_", "exp_N30_", "exp_N35_", "exp_N40_", "exp_N44_", "kNN_N30_", "kNN_N40_", "kNN_N50_", "spmv_N60_", "spmv_N65_", "spmv_N70_"};
const set<string> large = {"CG_N24_", "CG_N25_", "CG_N30_", "CG_N35_", "CG_N40_", "CG_N45_", "exp_N50_", "exp_N55_", "exp_N60_", "exp_N70_", "exp_N80_", "exp_N90_", "kNN_N45_", "kNN_N55_", "kNN_N60_", "kNN_N70_", "kNN_N75_", "kNN_N90_", "spmv_N120_", "spmv_N130_", "spmv_N150_"};


struct result
{
    string name;
    map<string, int> cost;
    map<string, int> supsteps;

    const set<string> inits = {"BSPg", "BSPgHCCS", "source3_weights_cluster", "source3_weights_clusterHCCS", "ILPInit", "ILPInitHCCS"};

    bool is(string s) { return name.find(s)!=string::npos;}

    bool has(string s) { return cost.find(s)!=cost.end();}

    bool in_set(const set<string>& S)
    {
        for(string s : S)
            if(is(s))
                return true;

        return false;
    }

    int getBestHeur()
    {
        string best = "None";
        int bestcost=-1;
        for(auto it = cost.begin(); it != cost.end(); ++it)
        {
            if(inits.find(it->first)==inits.end())
                continue;
            if(bestcost==-1 || it->second<bestcost)
            {
                best = it->first;
                bestcost = it->second;
            }
        }
        return bestcost;
    }

    int getRawInit()
    {
        string best = "None";
        int bestcost=-1;
        for(auto it = cost.begin(); it != cost.end(); ++it)
        {
            if(inits.find(it->first)==inits.end())
                continue;
            if(bestcost==-1 || it->second<bestcost)
            {
                best = it->first;
                bestcost = it->second;
            }
        }
        if(best.find("HCCS") == string::npos)
            return bestcost;
        string raw = best.substr(0, best.size()-4);
        if(cost.find(raw)!=cost.end())
            return cost[raw];
        cout<<"Error: raw version of algorithm not found for "<<name<<"."<<endl;
        return bestcost;
    }

    int getFinalCost()
    {
        if(cost.find("ILPCS")!=cost.end())
            return cost["ILPCS"];
        if(cost.find("ILPIter")!=cost.end())
            return cost["ILPIter"];
        if(cost.find("ILP")!=cost.end())
            return cost["ILP"];
        return getBestHeur();
    }

    int getILP()
    {
        if(cost.find("ILPIter")!=cost.end())
            return cost["ILPIter"];
        if(cost.find("ILP")!=cost.end())
            return cost["ILP"];
        return getBestHeur();
    }
    int getCompressed()
    {
        if(cost.find("ILPCS")!=cost.end())
            return cost["ILPCS"];
        if(cost.find("refined")!=cost.end())
            return cost["refined"];
        return -1;
    }
    int getHdagg()
    {
        if(cost.find("HDagg")!=cost.end())
            return cost["HDagg"];
        return -1;
    }
};

vector<result> readfile(string filename)
{
    ifstream infile(filename);
    if(!infile.is_open())
    {
        cout<<"Unable to find/open input file:"<<filename<<"\n";
        return vector<result>();
    }

    vector<result> results;
    string line;

    while(!infile.eof())
    {
        getline(infile, line);
        while(!infile.eof() && line.empty())
            getline(infile, line);
        if(infile.eof())
            break;

        result next;
        next.name = line;
        getline(infile, line);
        while(!infile.eof() && !line.empty())
        {
            stringstream stream(line);
            string temp;
            vector<string> values;
            while(getline(stream, temp, ','))
               values.push_back(temp);

            next.cost[values[0]] = stoi(values[1]);
            next.supsteps[values[0]] = stoi(values[2]);
            getline(infile, line);
            if(line.find(",") == string::npos)
                break;
        }
        results.push_back(next);
    }
    return results;
}

void addCompressed(vector<result>& results, string filename15, string filename30)
{
    vector<result> comp15 =readfile(filename15);
    vector<result> comp30 =readfile(filename30);

    // auxiliary structure to find names in comp15 and comp30 faster
    map<string, int> idx15;
    map<string, int> idx30;
    for(int i=0; i<comp15.size(); ++i)
        idx15[comp15[i].name]=i;
    for(int i=0; i<comp30.size(); ++i)
        idx30[comp30[i].name]=i;

    for(result& r : results)
    {
        int best = -1;
        if(idx15.find(r.name)!=idx15.end())
        {
            result r2 = comp15[idx15[r.name]];
            int cost = r2.getCompressed();
            if(best==-1 || cost<best)
                best=cost;
            r.cost["C15"] = cost;
        }
        if(idx30.find(r.name)!=idx30.end())
        {
            result r2 = comp30[idx30[r.name]];
            int cost = r2.getCompressed();
            if(best==-1 || cost<best)
                best=cost;
            r.cost["C30"] = cost;
        }
        if(best!=-1)
            r.cost["Copt"] = best;
    }
}

void addHDagg(vector<result>& results, string filenameHdagg)
{
    vector<result> hdagg = readfile(filenameHdagg);

    // auxiliary structure to find names in hdagg faster
    map<string, int> idx;
    for(int i=0; i<hdagg.size(); ++i)
        idx[hdagg[i].name]=i;

    for(result& r : results)
    {
        if(idx.find(r.name)!=idx.end())
        {
            result r2 = hdagg[idx[r.name]];
            r.cost["HDagg"] = r2.getHdagg();
        }
    }
}

int main()
{
    bool multi = true;
	// comment the following line if using the multilevel algorithm
    multi = false;

    // uncomment the line with the required file, or add a new one
    //vector<result> results =readfile("results_without_numa.csv");
    vector<result> results =readfile("results_with_numa.csv");
    //vector<result> results =readfile("results_huge.csv");
    //vector<result> results =readfile("results_huge_with_numa.csv");

    // also add the multilevel schedules to the results vector (coarsified to 15% and 30%), and HDagg
    if(multi)
        addCompressed(results, "results_multi15.csv", "results_multi30.csv");

    addHDagg(results, "results_hdagg.csv");

    // product for each algorithm
    double cilk=1, etf = 1, blest = 1, init = 1, heur = 1, ilp = 1, ilpcs = 1;
    // counter for each algorithm
    int cnt_cilk= 0, cnt_etf =0, cnt_blest = 0, cnt_init=0, cnt_heur=0, cnt_ilp=0, cnt_ilpcs=0;

    // product for multilevel variants
    double c15=1, c30 = 1, copt = 1;
    // counters for multilevel variants
    int cnt_c15= 0, cnt_c30 = 0, cnt_copt = 0;

    // product & counter for hdagg
    double hdagg = 1;
    int cnt_hdagg = 0;

    // Process cost ratios
    for(result r : results)
    {
        // uncomment one of these to normalize to HDagg, ETF or cilk
        //double base = (double) r.cost["HDagg"];
        //double base = (double) r.cost["ETF"];
        double base = (double) r.cost["cilk"];

        // uncomment some of the following to restrict the results to a given parameter or dataset
        /*if(!r.is("_numa4"))
            continue;*/
        /*if(!r.is("_p16"))
            continue;*/
        /*if(!r.is("_g1"))
            continue;*/
        /*if(!r.in_set(small))
            continue;*/

        if(multi && r.in_set(tiny))
            continue;

        //cout<<r.name<<endl;

        if(r.has("cilk"))
        {
            cilk *= (double) r.cost["cilk"] / base;
            ++cnt_cilk;
        }
        if(r.has("ETF"))
        {
            etf *= (double) r.cost["ETF"] / base;
            ++cnt_etf;
        }
        if(r.has("BL-EST"))
        {
            blest *= (double) r.cost["BL-EST"] / base;
            ++cnt_blest;
        }
        init *= (double) r.getRawInit() / base;
        heur *= (double) r.getBestHeur() / base;
        if(r.has("ILP") || r.has("ILPIter"))
        {
            ilp *= (double) r.getILP() / base;
            ++cnt_ilp;
        }
        ilpcs *= (double) r.getFinalCost() / base;

        ++cnt_init;
        ++cnt_heur;
        ++cnt_ilpcs;

        if(r.has("HDagg"))
        {
            hdagg *= (double) r.cost["HDagg"] / base;
            ++cnt_hdagg;
        }

        if(multi)
        {
            if(r.has("C15"))
            {
                c15 *= (double) r.cost["C15"] / base;
                ++cnt_c15;
            }
            if(r.has("C30"))
            {
                c30 *= (double) r.cost["C30"] / base;
                ++cnt_c30;
            }
            if(r.has("Copt"))
            {
                copt *= (double) r.cost["Copt"] / base;
                ++cnt_copt;
            }
        }

    }

    // Compute geometric mean
    cilk = pow(cilk,(double)1/cnt_cilk);
    etf = pow(etf,(double)1/cnt_etf);
    blest = pow(blest,(double)1/cnt_blest);
    init = pow(init,(double)1/cnt_init);
    heur = pow(heur,(double)1/cnt_heur);
    ilp = pow(ilp,(double)1/cnt_ilp);
    ilpcs = pow(ilpcs,(double)1/cnt_ilpcs);

    hdagg = pow(hdagg,(double)1/cnt_hdagg);

    if(multi)
    {
        c15 = pow(c15,(double)1/cnt_c15);
        c30 = pow(c30,(double)1/cnt_c30);
        copt = pow(copt,(double)1/cnt_copt);
    }

    // Print
    cout<<"Number of runs: "<<cnt_cilk<<endl;
    printf("Cilk: %.2lf \%\n", 100*cilk);
    printf("ETF: %.2lf \%\n", 100*etf);
    printf("BL-EST: %.2lf \%\n", 100*blest);
    printf("HDagg: %.2lf \%\n", 100*hdagg);
    printf("Init: %.2lf \%\n", 100*init);
    printf("HCCS: %.2lf \%\n", 100*heur);
    printf("ILP: %.2lf \%\n", 100*ilp);
    printf("Final: %.2lf \%\n\n", 100*ilpcs);

    if(multi)
    {
        printf("C15: %.2lf \%% (%d instances)\n", 100*c15, cnt_c15);
        printf("C30: %.2lf %% (%d instances)\n", 100*c30, cnt_c30);
        printf("Cx: %.2lf %% (%d instances)\n", 100*copt, cnt_copt);
    }

    return 0;
}
