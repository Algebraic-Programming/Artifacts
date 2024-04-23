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
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <set>
#include <map>
#include <climits>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

// List of initializaton methods available
vector<string> possibleModes{"random", "SJF", "cilk", "BSPg", "ETF", "BL-EST", "ETF-NUMA", "BL-EST-NUMA" };

// AUXILIARY FUNCTIONS

// unbiased random int generator
int randInt(int lim)
{
    int rnd = rand();
    while(rnd>=RAND_MAX-RAND_MAX%lim)
        rnd = rand();

    return rnd%lim;
}

struct intPair
{
    int a, b;

    intPair(){}
    intPair(int a1, int b1) : a(a1), b(b1) {}

    bool operator <(const intPair& other) const
    {
        return(a<other.a || (a==other.a && b<other.b));
    }
};

struct intTriple
{
    int a, b, c;
    intTriple(){}
    intTriple(int a1, int b1, int c1) : a(a1), b(b1), c(c1) {}
};


bool isDisjoint(vector<intPair>& intervals)
{
    if(intervals.empty())
        return true;

    sort(intervals.begin(), intervals.end());
    for(int i=0; i<intervals.size()-1; ++i)
        if(intervals[i].b>intervals[i+1].a)
            return false;

    return true;
}

int pow(int base, int exp)
{
    int val = 1;
    for(int i=0; i<exp; ++i)
        val *=base;
    return val;
}

// modify problem filename by adding substring at the right place
string editFilename(string filename, string toInsert)
{
    int pos = filename.find("_coarse");
    if(pos==string::npos)
        pos = filename.find("_instance");
    if(pos==string::npos)
        return toInsert + filename;

    return filename.substr(0, pos) + toInsert + filename.substr(pos, filename.length()-pos);
}

// PROBLEM-SPECIFIC

struct BSPproblem
{
    int p, g, L;
    vector<vector<int> > sendCost;
    int avgComm;
    BSPproblem(){}
    BSPproblem(int a, int b, int c) : p(a), g(b), L(c) { SetUniformCost(); }
    void SetUniformCost()
    {
        sendCost.clear();
        sendCost.resize(p, vector<int>(p, 1));
        for(int i=0; i<p; ++i)
            sendCost[i][i] = 0;
    }
    void SetExpCost(int base)
    {
        sendCost.clear();
        sendCost.resize(p, vector<int>(p, 0));
        int maxPos =1;
        for(; pow(2, maxPos+1)<=p-1; ++maxPos) ;
        for(int i=0; i<p; ++i)
            for(int j=i+1; j<p; ++j)
                for(int pos = maxPos; pos>=0; --pos)
                    if(((1<<pos)&i) != ((1<<pos)&j))
                    {
                        sendCost[i][j] = sendCost[j][i] = pow(base, pos);
                        break;
                    }
    }

	// write machine parameters to file
    void write(ofstream& outfile, bool NoNUMA = false) const
    {
        outfile<<p<<" "<<g<<" "<<L<<endl;

        if(!NoNUMA)
            for(int i=0; i<p; ++i)
                for(int j=0; j<p; ++j)
                    outfile<<i<<" "<<j<<" "<<sendCost[i][j]<<endl;

     }

	// compute average comm. coefficient between a pair of processors
     int computeCommAverage()
     {
         double avg = 0;
         for(int i=0; i<p; ++i)
            for(int j=0; j<p; ++j)
                avg += sendCost[i][j];
         avg = avg * (double) g / (double)p / (double)p;
         avgComm = round(avg);
         return avgComm;
     }
};

// main parameters for running simple_schedulers.cpp
struct ProgramParams
{
    bool NoNUMA = false;
    int TimeLimit = 600;
    bool NoHC = false;
    bool NoOutfile = false;
    bool OnlyHC = false;
    bool ContractMode = false;
    bool RefineMode = false;
    bool display = false;
    int newN = 0;
    int batchLength = 10;
    int HCsteps = 50;
    map<string, bool> runMode;
    ProgramParams()
    {
        for(string mode : possibleModes)
            runMode[mode]=false;
    }
};

// computational DAG
struct DAG
{
    int n;
	// in-neighbors and out-neighbors of each node
    vector< vector <int> > In, Out;

	// work and communication weight of each node
    vector<int> workW, commW;

    void Resize(int N)
    {
        n=N;
        In.clear();
        In.resize(n);
        Out.clear();
        Out.resize(n);
        workW.clear();
        workW.resize(n);
        commW.clear();
        commW.resize(n);
    }

    void addEdge(int v1, int v2, bool noPrint = false)
    {
        if(v1>=v2)
            cout<<"DAG edge addition error."<<endl;

        if(v2>=n)
            cout<<"Error: node index out of range."<<endl;

        In[v2].push_back(v1);
        Out[v1].push_back(v2);
    }

	// compute topological ordering
    vector<int> GetTopOrder() const
    {
        vector<int> predecessors(n, 0);
        list<int> next;

        vector<int> TopOrder;

        // Find source nodes
        for(int i=0; i<n; ++i)
            if(In[i].empty())
                next.push_back(i);

        // Execute BFS
        while(!next.empty())
        {
            int node = next.front();
            next.pop_front();
            TopOrder.push_back(node);

            for(int i=0; i<Out[node].size(); ++i)
            {
                int current = Out[node][i];
                ++predecessors[current];
                if(predecessors[current]==In[current].size())
                    next.push_back(current);
            }
        }

        if(TopOrder.size()!=n)
            cout<<"Error during topological ordering!"<<endl;

        return TopOrder;
    }

	// read DAG from file in hyperDAG format
    bool read(ifstream& infile)
    {
        string line;
        getline(infile, line);
        while(!infile.eof() && line.at(0)=='%')
           getline(infile, line);

        int hEdges, pins, N;
        sscanf(line.c_str(), "%d %d %d", &hEdges, &N, &pins);

        if(N<=0 || hEdges<=0 || pins<=0)
        {
            cout<<"Incorrect input file format (number of nodes/hyperedges/pins is not positive).\n";
            return false;
        }

        Resize(N);
        vector<int> edgeSource(hEdges, -1);
        // read edges
        for(int i=0; i<pins; ++i)
        {
            if(infile.eof())
            {
                cout<<"Incorrect input file format (file terminated too early).\n";
                return false;
            }
            getline(infile, line);
            while(!infile.eof() && line.at(0)=='%')
                getline(infile, line);

            int hEdge, node;
            sscanf(line.c_str(), "%d %d", &hEdge, &node);

            if(hEdge<0 || node<0 || hEdge>=hEdges || node>=N)
            {
                cout<<"Incorrect input file format (index out of range).\n";
                return false;
            }

            if(edgeSource[hEdge]==-1)
                edgeSource[hEdge]=node;
            else
                addEdge(edgeSource[hEdge], node);
        }

        ReOrderDAGVectors();

        for(int i=0; i<N; ++i)
        {
            if(infile.eof())
            {
                cout<<"Incorrect input file format (file terminated too early).\n";
                return false;
            }

            getline(infile, line);
            while(!infile.eof() && line.at(0)=='%')
                getline(infile, line);

            int node, work, comm;
            sscanf(line.c_str(), "%d %d %d", &node, &work, &comm);

            if(node<0 || work<0 || comm<0 || node>=N)
            {
                cout<<"Incorrect input file format (index out of range, our weight below 0).\n";
                return false;
            }

            workW[node]=work;
            commW[node]=comm;
        }

        return true;
    }

    bool read(string filename)
    {
        ifstream infile(filename);
        if(!infile.is_open())
        {
            cout<<"Unable to find/open input schedule file.\n";
            return false;
        }

        read(infile);
        infile.close();
        return true;
    }

	// write DAG to file in hyperDAG format
    void write(ofstream& outfile)
    {
        int sinks=0, pins=0;
        for(int i=0; i<n; ++i)
            if(Out[i].size()>0)
                pins+=1+Out[i].size();
            else
                ++sinks;

        outfile << n - sinks << " "<< n << " " << pins <<"\n";

        int edgeIndex = 0;
        for(int i=0; i<n; ++i)
            if(Out[i].size()>0)
            {
                outfile << edgeIndex << " "<< i <<"\n";
                for(int j=0; j<Out[i].size(); ++j)
                    outfile << edgeIndex << " "<< Out[i][j] <<"\n";

                ++edgeIndex;
            }

        for(int i=0; i<n; ++i)
            outfile << i << " "<< workW[i] << " " << commW[i] <<"\n";

     }

	// ensure that the list of in- and out-neighbors for nodes follows a topological ordering
    void ReOrderDAGVectors()
    {
        vector<int> topOrder = GetTopOrder();
        vector< vector <int> > newIn, newOut;
        newIn.resize(n);
        newOut.resize(n);
        for(int i=0; i<n; ++i)
        {
            int node = topOrder[i];
            for(int j : In[node])
                newOut[j].push_back(node);
            for(int j : Out[node])
                newIn[j].push_back(node);
        }
        In = newIn;
        Out = newOut;
    }

	// get topological order for a subset of (valid) nodes
    vector<int> GetCleanedTopOrder(const vector<bool>& valid) const
    {
        vector<int> TopOrder = GetTopOrder(), cleanedOrder;
        for(int node : TopOrder)
            if(valid[node])
                cleanedOrder.push_back(node);

        return cleanedOrder;
    }

	// main functions for multilevel approach
    DAG Coarsify(int newN, string outfile, int randomSeed) const;
    DAG Contract(const vector<intPair>& contractionSteps) const;

	// auxiliary functions for multilevel approach
    map<intPair, bool> GetContractableEdges() const;
    bool isContractable(int source, int target, const vector<int>& topOrder, const vector<bool>& valid) const;
    void updateDistantEdgeContractibility(int source, int target, map<intPair, bool>& contractable) const;
    int getLongestPath(const set<int>& nodes) const;

};

// TOOLS FOR MULTILEVEL SCHEDULING

// check if a DAG edge is contractable without creating directed cycles
// (i.e. if there is directed path from source node to target node besides the edge)
// (further inputs: the topological order index of each node, and an indicator of which nodes are still valid)
bool DAG::isContractable(int source, int target, const vector<int>& topOrderPos, const vector<bool>& valid) const
{
    list<int> Queue;
    set<int> visited;
    for(int succ: Out[source])
        if(valid[succ] && topOrderPos[succ]<topOrderPos[target])
        {
            Queue.push_back(succ);
            visited.insert(succ);
        }

    while(!Queue.empty())
    {
        int node = Queue.front();
        Queue.pop_front();
        for(int succ: Out[node])
        {
            if(succ == target)
                return false;

            if(valid[succ] && topOrderPos[succ]<topOrderPos[target] && visited.count(succ)==0)
            {
                Queue.push_back(succ);
                visited.insert(succ);
            }
        }
    }
    return true;
}

// determine contractability for each edge of the DAG
map<intPair, bool> DAG::GetContractableEdges() const
{
    vector<int> topOrder = GetTopOrder();
    vector<int> topOrderPos(n);
    for(int i=0; i<topOrder.size(); ++i)
            topOrderPos[topOrder[i]] = i;

    vector<bool> valid(n, true);
    map<intPair, bool> contractable;

    for(int i=0; i<n; ++i)
        for(int succ : Out[i])
            contractable[intPair(i, succ)] = isContractable(i, succ, topOrderPos, valid);

    return contractable;
}

// after the contraction of an edge, update the contractability of other edges in the DAG
void DAG::updateDistantEdgeContractibility(int source, int target, map<intPair, bool>& contractable) const
{
    vector<bool> descendant(n, false);
    set<int> ancestors;
    list<int> Queue;

    for(int succ: Out[source])
        if(succ != target)
        {
            Queue.push_back(succ);
            descendant[succ]=true;
        }
    while(!Queue.empty())
    {
        int node = Queue.front();
        Queue.pop_front();
        for(int succ: Out[node])
            if(!descendant[succ])
            {
                Queue.push_back(succ);
                descendant[succ] = true;
            }
    }

    for(int pred: In[target])
        if(pred != source)
        {
            Queue.push_back(pred);
            ancestors.insert(pred);
        }
    while(!Queue.empty())
    {
        int node = Queue.front();
        Queue.pop_front();
        for(int pred: In[node])
            if(ancestors.count(pred)==0)
            {
                Queue.push_back(pred);
                ancestors.insert(pred);
            }
    }

    for(int node : ancestors)
        for(int succ : Out[node])
            if(descendant[succ])
                contractable[intPair(node, succ)] = false;
}

// get longest path in an induced subgraph
int DAG::getLongestPath(const set<int>& nodes) const
{
    list<int> Q;
    map<int, int> dist, inDegree, visited;

    // Find source nodes
    for(int node : nodes)
    {
        int indeg;
        for(int pred : In[node])
            if(nodes.count(pred)==1)
                ++indeg;

        if(indeg==0)
        {
            Q.push_back(node);
            dist[node] = 0;
        }
        inDegree[node]=indeg;
        visited[node]=0;
    }

    // Execute BFS
    while(!Q.empty())
    {
        int node = Q.front();
        Q.pop_front();

        for(int succ : Out[node])
        {
            if(nodes.count(succ)==0)
                continue;

            ++visited[succ];
            if(visited[succ]==inDegree[succ])
            {
                Q.push_back(succ);
                dist[succ] = dist[node] + 1;
            }
        }
    }

    int mx = 0;
    for(int node : nodes)
        mx = max(mx, dist[node]);

    return mx;
}

// auxiliary structure to store contractable edges
struct contractionEdge
{
    intPair edge;
    int nodeW;
    int edgeW;

    contractionEdge(int from, int to, int Wnode, int Wedge) : edge(from, to), nodeW(Wnode), edgeW(Wedge) {}

    bool operator <(const contractionEdge& other) const
    {
        return(nodeW<other.nodeW || (nodeW==other.nodeW && edgeW<other.edgeW));
    }
};


// compute the image (new index) of each original node after a list of contraction steps
vector<int> GetFinalImage(const DAG& G, const vector<intPair>& contractionSteps)
{
    vector<int> target(G.n), pointsTo(G.n, -1);

    for(const intPair& step : contractionSteps)
        pointsTo[step.b] = step.a;


    for(int i=0; i<G.n; ++i)
    {
        target[i]=i;
        while(pointsTo[target[i]]!=-1)
            target[i] = pointsTo[target[i]];
    }

    vector<bool> valid(G.n, false);
    for(int i=0; i<G.n; ++i)
        valid[target[i]]=true;

    DAG G2;
    G2.Resize(G.n);
    set<intPair> edges;
    for(int i=0; i<G.n; ++i)
        for(int succ : G.Out[i])
        {
            intPair edge(target[i], target[succ]);
            if(edges.find(edge)!=edges.end() || edge.a==edge.b)
                continue;
            G2.Out[edge.a].push_back(edge.b);
            G2.In[edge.b].push_back(edge.a);
            edges.insert(edge);
        }

    vector<int> topOrder = G2.GetCleanedTopOrder(valid);

    vector<int> newIdx(G.n);
    for(int i=0; i<topOrder.size(); ++i)
    {
        if(target[topOrder[i]]!=topOrder[i])
            cout<<"ERROR: Invalid contraction output."<<endl;

        newIdx[topOrder[i]] = i;
    }

    for(int i=0; i<G.n; ++i)
        target[i] = newIdx[target[i]];

    return target;
}


// main function: coarsify input DAG into newN nodes
DAG DAG::Coarsify(int newN, string outfilename, int randomSeed = -1) const
{
    // Init edge weights
    map<intPair, int> edges;
    for(int i=0; i<n; ++i)
        for(int j : Out[i])
           edges[intPair(i, j)] = commW[i];

	// Init topological order
    vector<int> topOrder = GetTopOrder();
    vector<int> topOrderPos(n);
    for(int i=0; i<topOrder.size(); ++i)
        topOrderPos[topOrder[i]] = i;

    vector<bool> validNode(n, true);
    vector<intPair> contractionHistory, topOrderHistory;
    map<intPair, bool> contractable = GetContractableEdges();

    map<intPair, int> longestPath;
    for(auto it = contractable.begin(); it != contractable.end(); ++it)
        longestPath[it->first] = 1;

	// list of original node indices contained in each contracted node
    vector<set<int> > contains(n);
    for(int i=0; i<n; ++i)
        contains[i].insert(i);


    DAG G = *this;

    bool random = (randomSeed != -1);
    if(random)
        srand(randomSeed);

    for(int NrOfNodes = n; NrOfNodes > newN; --NrOfNodes)
    {
        // Single contraction step

        //Collect contraction edge candidates
        vector<contractionEdge> candidates;
        for(int i=0; i<n; ++i)
        {
            if(!validNode[i] || G.Out[i].empty())
                continue;

            for(int succ : G.Out[i])
                if(validNode[succ] && contractable[intPair(i, succ)])
                    candidates.push_back(contractionEdge(i, succ, contains[i].size()+contains[succ].size(), edges[intPair(i, succ)]));
        }

        if(candidates.empty())
        {
            cout<<"ERROR: no edge to contract. "<<NrOfNodes<<endl;
            break;
        }

        //Select edge to contract
        sort(candidates.begin(), candidates.end());
        int limit = (candidates.size()+2)/3;
        int limitWeight = candidates[limit].nodeW;
        while(limit<candidates.size()-1 && candidates[limit+1].nodeW == limitWeight)
            ++limit;

        // an edge case
        if(candidates.size()==1)
            limit = 0;

        contractionEdge next = candidates[randInt(candidates.size())];
        if(!random)
        {
            int best = 0;
            for(int i=1; i<=limit; ++i)
                if(candidates[i].edgeW > candidates[best].edgeW)
                    best = i;

            next = candidates[best];
        }

        //Find far-away edges that become uncontractable now
        G.updateDistantEdgeContractibility(next.edge.a, next.edge.b, contractable);

        //Contract edge
        G.workW[next.edge.a] += G.workW[next.edge.b];
        G.commW[next.edge.a] += G.commW[next.edge.b];
        validNode[next.edge.b] = false;
        contractionHistory.push_back(intPair(next.edge.a, next.edge.b));
        for(int pred : G.In[next.edge.b])
        {
            if(pred == next.edge.a)
                continue;

            int idx = -1;
            for(int j=0; j< G.In[next.edge.a].size(); ++j)
                if(G.In[next.edge.a][j]==pred)
                {
                    idx = j;
                    break;
                }

            if(idx >= 0) // Combine edges
            {
                edges[intPair(pred, next.edge.a)] = 0;
                for(int node : contains[pred])
                    for(int succ : G.Out[node])
                        if(succ == next.edge.a || succ == next.edge.b)
                            edges[intPair(pred, next.edge.a)] += commW[node];
            }
            else // Add incoming edge
            {
                edges[intPair(pred, next.edge.a)] = edges[intPair(pred, next.edge.b)];
                G.Out[pred].push_back(next.edge.a);
                G.In[next.edge.a].push_back(pred);
            }
            for(auto it = G.Out[pred].begin(); it != G.Out[pred].end(); ++it)
                if(*it==next.edge.b)
                {
                    G.Out[pred].erase(it);
                    break;
                }
        }
        for(int succ : G.Out[next.edge.b])
        {
            int idx = -1;
            for(int j=0; j< G.Out[next.edge.a].size(); ++j)
                if(G.Out[next.edge.a][j]==succ)
                {
                    idx = j;
                    break;
                }

            if(idx >= 0) // Combine edges
            {
                edges[intPair(next.edge.a, succ)] += edges[intPair(next.edge.b, succ)];
            }
            else // Add outgoing edge
            {
                edges[intPair(next.edge.a, succ)] = edges[intPair(next.edge.b, succ)];
                G.In[succ].push_back(next.edge.a);
                G.Out[next.edge.a].push_back(succ);
            }
            for(auto it = G.In[succ].begin(); it != G.In[succ].end(); ++it)
                if(*it==next.edge.b)
                {
                    G.In[succ].erase(it);
                    break;
                }
        }
        for(auto it = G.Out[next.edge.a].begin(); it != G.Out[next.edge.a].end(); ++it)
            if(*it==next.edge.b)
            {
                G.Out[next.edge.a].erase(it);
                break;
            }

        G.In[next.edge.b].clear();
        G.Out[next.edge.b].clear();

        for(int node : contains[next.edge.b])
            contains[next.edge.a].insert(node);

        contains[next.edge.b].clear();

        //Update topological order
        topOrder = G.GetCleanedTopOrder(validNode);
        for(int i=0; i<topOrder.size(); ++i)
            topOrderPos[topOrder[i]] = i;

        //Update contractable edges
        for(int pred : G.In[next.edge.a])
        {
            contractable[intPair(pred, next.edge.a)] = G.isContractable(pred, next.edge.a, topOrderPos, validNode);
            if(contractable[intPair(pred, next.edge.a)])
            {
                set<int> merged = contains[pred];
                for(int node : contains[next.edge.a])
                    merged.insert(node);

                longestPath[intPair(pred, next.edge.a)] = getLongestPath(merged);
            }
        }

        for(int succ : G.Out[next.edge.a])
        {
            contractable[intPair(next.edge.a, succ)] = G.isContractable(next.edge.a, succ, topOrderPos, validNode);
            if(contractable[intPair(next.edge.a, succ)])
            {
                set<int> merged = contains[succ];
                for(int node : contains[next.edge.a])
                    merged.insert(node);

                longestPath[intPair(next.edge.a, succ)] = getLongestPath(merged);
            }
        }


    }

    //Print contraction steps to file
    ofstream outfile(outfilename);
    if(outfile.is_open())
    {
        outfile<<n<<" "<<newN<<endl;
        for(intPair entry : contractionHistory)
            outfile<<entry.a<<" "<<entry.b<<endl;

        vector<int> image = GetFinalImage(*this, contractionHistory);
        for(int i=0; i<n; ++i)
            outfile<<i<<" "<<image[i]<<endl;

        outfile.close();
    }
    else
        cout<<"ERROR: Unable to write/open contraction output log file for DAG.\n";

    //Return contracted DAG
    return Contract(contractionHistory);

}

// Given the selected contraction steps, return the contracted DAG
DAG DAG::Contract(const vector<intPair>& contractionSteps) const
{
    vector<int> target = GetFinalImage(*this, contractionSteps);

    DAG G;
    G.Resize(n-contractionSteps.size());
    for(int i=0; i<G.n; ++i)
    {
        G.workW[i] = 0;
        G.commW[i] = 0;
    }

    for(int i=0; i<n; ++i)
    {
        G.workW[target[i]] += workW[i];
        G.commW[target[i]] += commW[i];
    }

    set<intPair> edges;
    for(int i=0; i<n; ++i)
        for(int succ : Out[i])
        {
            intPair edge(target[i], target[succ]);
            if(edges.find(edge)!=edges.end() || edge.a==edge.b)
                continue;
            G.addEdge(edge.a, edge.b);
            edges.insert(edge);
        }

    for(int i=0; i<G.n; ++i)
    {
        sort(G.In[i].begin(), G.In[i].end());
        sort(G.Out[i].begin(), G.Out[i].end());
    }

    return G;
}

// Coarsify DAG, and write result to file
bool CoarsifyAndWrite(const DAG& G, const BSPproblem& params, int newN, string filename, string contractFile, bool NoNUMA = false)
{
    DAG coarse = G.Coarsify(newN, contractFile, -1);

    ofstream outfile(filename);
    if(!outfile.is_open())
    {
        cout<<"Unable to write/open problem file for coarsified DAG.\n";
        return false;
    }
    coarse.write(outfile);
    params.write(outfile, NoNUMA);
    outfile.close();
    return true;
}

// Process a file with list of contraction steps
bool ReadContractionFile(string filename, vector<intPair>& contractionSteps)
{
    ifstream infile(filename);
    if(!infile.is_open())
    {
        cout<<"Unable to find/open input DAG-contraction file.\n";
        return false;
    }

    string line;
    getline(infile, line);
    while(!infile.eof() && line.at(0)=='%')
       getline(infile, line);

    int oldN, newN;
    sscanf(line.c_str(), "%d %d", &oldN, &newN);

    contractionSteps.clear();

    for(int i=0; i<oldN-newN; ++i)
    {
        if(infile.eof())
        {
            cout<<"Incorrect contraction file format (file terminated too early).\n";
            return false;
        }
        getline(infile, line);
        while(!infile.eof() && line.at(0)=='%')
            getline(infile, line);

        int target, source;
        sscanf(line.c_str(), "%d %d", &target, &source);

        if(target<0 || source<0 || target>=oldN || source>=oldN)
        {
            cout<<"Incorrect contraction file format (index out of range).\n";
            return false;
        }

        contractionSteps.push_back(intPair(target, source));
    }

    infile.close();
    return true;
}


// BSP PROBLEM PARAMETER HANDLING

// Reading machine parameters from file
bool readProblemParams(ifstream& infile, BSPproblem& params, bool NoNUMA=false)
{
    string line;
    getline(infile, line);
    while(!infile.eof() && line.at(0)=='%')
        getline(infile, line);

    sscanf(line.c_str(), "%d %d %d", &params.p, &params.g, &params.L);

    params.SetUniformCost();

    if(!NoNUMA)
    {
        for(int i=0; i<params.p*params.p; ++i)
        {
            if(infile.eof())
            {
                cout<<"Incorrect input file format (file terminated too early).\n";
                return false;
            }
            getline(infile, line);
            while(!infile.eof() && line.at(0)=='%')
                getline(infile, line);

            int fromProc, toProc, value;
            sscanf(line.c_str(), "%d %d %d", &fromProc, &toProc, &value);

            if(fromProc<0 || toProc<0 || fromProc>=params.p || toProc>=params.p || value<0)
            {
                cout<<"Incorrect input file format (index out of range or negative NUMA value).\n";
                return false;
            }
            if(fromProc==toProc && value!=0)
            {
                cout<<"Incorrect input file format (main diagonal of NUMA cost matrix must be 0).\n";
                return false;
            }
            params.sendCost[fromProc][toProc]=value;
        }
    }
    params.computeCommAverage();

    return true;
}

bool readProblem(string filename, DAG& G, BSPproblem& params, bool NoNUMA=true)
{
    ifstream infile(filename);
    if(!infile.is_open())
    {
        cout<<"Unable to find/open input problem file.\n";
        return false;
    }

    G.read(infile);
    readProblemParams(infile, params, NoNUMA);

    infile.close();
    return true;
}


// MAIN STRUCTURE: SCHEDULE
// supports both BSP and classical schedules

struct schedule
{
	// scheduling problem
    DAG G;
    BSPproblem params;

	// main ingredients of a schedule
	// (note: time is for time points in classical schedule; supstep is for supersteps in BSP schedule)
    vector<int> proc, time, supstep;
    vector<vector<list<int> > > supsteplists;
	int cost;

    //comm schedule
    vector<vector<int > > commSchedule;

    // aux data for greedy scheduling algorithms
    vector<vector<int> > procInHyperedge;
    vector<list<int> > procQueue;
    vector<list<int> > greedyProcLists;

    // aux data for hill climbing
    vector<vector<vector<bool> > > canMove;
    vector<list<intPair> > moveOptions;
    vector<vector<vector<list<intPair>::iterator> > > movePointer;
    vector<vector<map<int, int> > > succSteps;
    vector<vector<int> > workCost, sent, received, commCost;
    vector<set< intPair > > workCostList, commCostList;
    vector<vector<set<intPair>::iterator> > workCostPointer, commCostPointer;
    vector<list<int>::iterator > supStepListPointer;
    bool HCwithLatency = true;

	// aux data for comm schedule hill climbing
    vector<vector<intPair > > commBounds;
    vector<vector<list<intPair> > > commSchedSendLists;
    vector<vector<list<intPair>::iterator > > commSchedSendListPointer;
    vector<vector<list<intPair> > > commSchedRecLists;
    vector<vector<list<intPair>::iterator > > commSchedRecListPointer;

    enum Direction {
        EARLIER = 0,
        AT,
        LATER
    };
    static const int NumDirections = 3;

    struct stepAuxData
    {
        int newCost;
        map<intPair, int> sentChange, recChange;
        bool canShrink = false;
    };

	// INGREDIENTS OF HILL CLIMBING

	// Initialize data structures (based on current schedule)
    void InitHillClimbing()
    {
        int N = G.n;
        int M = supsteplists.size();
        cost = GetBSPCost();

        //Movement options
        canMove.clear();
        canMove.resize(NumDirections, vector<vector<bool> >(N, vector<bool>(params.p, false)));
        moveOptions.clear();
        moveOptions.resize(NumDirections);
        movePointer.clear();
        movePointer.resize(NumDirections, vector<vector<list<intPair>::iterator> >(N, vector<list<intPair>::iterator>(params.p)));

        //Value use lists
        succSteps.clear();
        succSteps.resize(N, vector<map<int, int> >(params.p));
        for(int i=0; i<N; ++i)
            for(int succ : G.Out[i])
            {
                if(succSteps[i][proc[succ]].find(supstep[succ])==succSteps[i][proc[succ]].end())
                    succSteps[i][proc[succ]].insert({supstep[succ], 1});
                else
                    succSteps[i][proc[succ]].at(supstep[succ]) += 1;
            }

        //Cost data
        workCost.clear();
        workCost.resize(M, vector<int>(params.p, 0));
        sent.clear();
        sent.resize(M-1, vector<int>(params.p, 0));
        received.clear();
        received.resize(M-1, vector<int>(params.p, 0));
        commCost.clear();
        commCost.resize(M-1, vector<int>(params.p));

        workCostList.clear();
        workCostList.resize(M);
        commCostList.clear();
        commCostList.resize(M-1);
        workCostPointer.clear();
        workCostPointer.resize(M, vector<set<intPair>::iterator>(params.p));
        commCostPointer.clear();
        commCostPointer.resize(M-1, vector<set<intPair>::iterator>(params.p));

        //Supstep list pointers
        supStepListPointer.clear();
        supStepListPointer.resize(N);
        for(int i=0; i<M; ++i)
            for(int j=0; j<params.p; ++j)
                for(auto it = supsteplists[i][j].begin(); it!=supsteplists[i][j].end(); ++it)
                    supStepListPointer[*it] = it;

        //Compute movement options
        for(int i=0; i<N; ++i)
            updateNodeMoves(i);

        //Compute cost data
        cost=0;
        for(int i=0; i<M; ++i)
        {
            for(int j=0; j<params.p; ++j)
            {
                for(int node : supsteplists[i][j])
                    workCost[i][j] += G.workW[node];

                intPair entry(workCost[i][j], j);
                workCostPointer[i][j] = workCostList[i].insert(entry).first;
            }
            cost+=(--workCostList[i].end())->a;
        }

        vector<vector <bool> > present(N, vector<bool>(params.p, false));
        for(int i=0; i<M-1; ++i)
        {
            for(int j=0; j<params.p; ++j)
                for(int node : supsteplists[i+1][j])
                    for(int pred : G.In[node])
                        if(proc[node]!=proc[pred] && !present[pred][proc[node]])
                        {
                            present[pred][proc[node]] = true;
                            sent[i][proc[pred]] += G.commW[pred] * params.sendCost[proc[pred]][proc[node]];
                            received[i][proc[node]] += G.commW[pred] * params.sendCost[proc[pred]][proc[node]];
                        }

            for(int j=0; j<params.p; ++j)
            {
                commCost[i][j]=max(sent[i][j], received[i][j]);
                intPair entry(commCost[i][j], j);
                commCostPointer[i][j] = commCostList[i].insert(entry).first;
            }
            cost+=params.g * commCostList[i].rbegin()->a + params.L;
        }
    }

	// Functions to ompute and update the list of possible moves
    void updateNodeMovesEarlier(int node)
    {
        if(supstep[node]==0)
            return;

        set<int> predProc;
        for(int pred : G.In[node])
        {
            if(supstep[pred]==supstep[node])
                return;
            if(supstep[pred]==supstep[node]-1)
                predProc.insert(proc[pred]);
        }

        if(predProc.size()>1)
            return;

        if(predProc.size()==1)
            addMoveOption(node, *predProc.begin(), EARLIER);
        else for(int j=0; j<params.p; ++j)
            addMoveOption(node, j, EARLIER);
    }

    void updateNodeMovesAt(int node)
    {
        for(int pred : G.In[node])
            if(supstep[pred]==supstep[node])
                return;

        for(int succ : G.Out[node])
            if(supstep[succ]==supstep[node])
                return;

        for(int j=0; j<params.p; ++j)
            if(j!=proc[node])
                addMoveOption(node, j, AT);
    }

    void updateNodeMovesLater(int node)
    {
        if(supstep[node]==supsteplists.size()-1)
            return;

        set<int> succProc;
        for(int succ : G.Out[node])
        {
            if(supstep[succ]==supstep[node])
                return;
            if(supstep[succ]==supstep[node]+1)
                succProc.insert(proc[succ]);
        }

        if(succProc.size()>1)
            return;

        if(succProc.size()==1)
            addMoveOption(node, *succProc.begin(), LATER);
        else for(int j=0; j<params.p; ++j)
            addMoveOption(node, j, LATER);
    }

    void updateNodeMoves(int node)
    {
        eraseMoveOptions(node);
        updateNodeMovesEarlier(node);
        updateNodeMovesAt(node);
        updateNodeMovesLater(node);
    }

    void addMoveOption(int node, int p, Direction dir)
    {
        if(!canMove[dir][node][p])
        {
            canMove[dir][node][p] = true;
            moveOptions[dir].push_back(intPair(node, p));
            movePointer[dir][node][p] = --moveOptions[dir].end();
        }
    }

    void eraseMoveOptions(int node)
    {
        for(int j=0; j<params.p; ++j)
            for(int dir = 0; dir<NumDirections; ++dir)
                if(canMove[dir][node][j])
                {
                    canMove[dir][node][j] = false;
                    moveOptions[dir].erase(movePointer[dir][node][j]);
                }
    }

	// Compute the cost change incurred by a potential move
    int moveCostChange(int node, int p, int where, stepAuxData& changing)
    {
        int step = supstep[node];
        int oldProc = proc[node];
        int change = 0;

        // Work cost change
        auto itBest = --workCostList[step].end();
        int maxAfterRemoval = itBest->a;
        if(itBest->b == oldProc)
        {
            auto itNext = itBest;
            --itNext;
            maxAfterRemoval = max(itBest->a-G.workW[node], itNext->a);
            change -= itBest->a - maxAfterRemoval;
        }

        int maxBeforeAddition = (where==0) ? maxAfterRemoval : workCostList[step+where].rbegin()->a;
        if(workCost[step+where][p]+G.workW[node] > maxBeforeAddition)
            change += workCost[step+where][p]+G.workW[node] - maxBeforeAddition;

        // Comm cost change
        list<intTriple> sentInc, recInc;
        //  -outputs
        if(p!=oldProc)
        {
            for(int j=0; j<params.p; ++j)
            {
                if(succSteps[node][j].empty())
                    continue;

                int affectedStep = succSteps[node][j].begin()->first - 1;
                if(j==p)
                {
                    sentInc.push_back(intTriple(affectedStep, oldProc, -G.commW[node]*params.sendCost[oldProc][j]));
                    recInc.push_back(intTriple(affectedStep, p, -G.commW[node]*params.sendCost[oldProc][j]));
                }
                else if(j==oldProc)
                {
                    recInc.push_back(intTriple(affectedStep, oldProc, G.commW[node]*params.sendCost[p][j]));
                    sentInc.push_back(intTriple(affectedStep, p, G.commW[node]*params.sendCost[p][j]));
                }
                else
                {
                    sentInc.push_back(intTriple(affectedStep, oldProc, -G.commW[node]*params.sendCost[oldProc][j]));
                    recInc.push_back(intTriple(affectedStep, j, -G.commW[node]*params.sendCost[oldProc][j]));
                    sentInc.push_back(intTriple(affectedStep, p, G.commW[node]*params.sendCost[p][j]));
                    recInc.push_back(intTriple(affectedStep, j, G.commW[node]*params.sendCost[p][j]));
                }
            }
        }

        //  -inputs
        if(p==oldProc)
            for(int pred : G.In[node])
            {
                if(proc[pred]==p)
                    continue;

                auto firstUse = *succSteps[pred][p].begin();
                bool skip = firstUse.first<step || (firstUse.first == step && where>=0 && firstUse.second>1);
                if(!skip)
                {
                    sentInc.push_back(intTriple(step-1, proc[pred], -G.commW[pred]*params.sendCost[proc[pred]][p]));
                    recInc.push_back(intTriple(step-1, p, -G.commW[pred]*params.sendCost[proc[pred]][p]));
                    sentInc.push_back(intTriple(step+where-1, proc[pred], G.commW[pred]*params.sendCost[proc[pred]][p]));
                    recInc.push_back(intTriple(step+where-1, p, G.commW[pred]*params.sendCost[proc[pred]][p]));
                }
            }
        else
            for(int pred : G.In[node])
            {
                // Comm. cost of sending pred to oldProc
                auto firstUse = succSteps[pred][oldProc].begin();
                bool skip = (proc[pred]==oldProc) || firstUse->first<step || (firstUse->first== step && firstUse->second>1);
                if(!skip)
                {
                    sentInc.push_back(intTriple(step-1, proc[pred], -G.commW[pred]*params.sendCost[proc[pred]][oldProc]));
                    recInc.push_back(intTriple(step-1, oldProc, -G.commW[pred]*params.sendCost[proc[pred]][oldProc]));
                    ++firstUse;
                    if(firstUse!=succSteps[pred][oldProc].end())
                    {
                        int nextStep = firstUse->first;
                        sentInc.push_back(intTriple(nextStep-1, proc[pred], G.commW[pred]*params.sendCost[proc[pred]][oldProc]));
                        recInc.push_back(intTriple(nextStep-1, oldProc, G.commW[pred]*params.sendCost[proc[pred]][oldProc]));
                    }
                }

                // Comm. cost of sending pred to p
                firstUse = succSteps[pred][p].begin();
                skip = (proc[pred]==p) || ((firstUse != succSteps[pred][p].end()) && (firstUse->first <= step+where));
                if(!skip)
                {
                    sentInc.push_back(intTriple(step+where-1, proc[pred], G.commW[pred]*params.sendCost[proc[pred]][p]));
                    recInc.push_back(intTriple(step+where-1, p, G.commW[pred]*params.sendCost[proc[pred]][p]));
                    if(firstUse != succSteps[pred][p].end())
                    {
                        sentInc.push_back(intTriple(firstUse->first-1, proc[pred], -G.commW[pred]*params.sendCost[proc[pred]][p]));
                        recInc.push_back(intTriple(firstUse->first-1, p, -G.commW[pred]*params.sendCost[proc[pred]][p]));
                    }
                }

            }

        //  -process changes
        changing.sentChange.clear();
        changing.recChange.clear();
        set<int> affectedSteps;
        for(auto entry : sentInc)
        {
            affectedSteps.insert(entry.a);
            auto itr = changing.sentChange.find(intPair(entry.a, entry.b));
            if(itr==changing.sentChange.end())
                changing.sentChange.insert({intPair(entry.a, entry.b), entry.c});
            else
                itr->second += entry.c;
        }
        for(auto entry : recInc)
        {
            affectedSteps.insert(entry.a);
            auto itr = changing.recChange.find(intPair(entry.a, entry.b));
            if(itr==changing.recChange.end())
                changing.recChange.insert({intPair(entry.a, entry.b), entry.c});
            else
                itr->second += entry.c;
        }

        auto itrSent = changing.sentChange.begin(), itrRec = changing.recChange.begin();
        for(int sstep : affectedSteps)
        {
            int newMax=0;
            for(int j=0; j<params.p; ++j)
            {
                int diff = (itrSent!=changing.sentChange.end() && itrSent->first.a==sstep && itrSent->first.b==j) ? (itrSent++)->second : 0;
                if(sent[sstep][j]+diff>newMax)
                    newMax = sent[sstep][j]+diff;
                diff = (itrRec != changing.recChange.end() && itrRec->first.a==sstep && itrRec->first.b==j) ? (itrRec++)->second : 0;
                if(received[sstep][j]+diff>newMax)
                    newMax = received[sstep][j]+diff;
            }
            change += params.g * (newMax - commCostList[sstep].rbegin()->a);

            if(HCwithLatency)
            {
                if(newMax>0 && commCostList[sstep].rbegin()->a==0)
                {
                    change += params.L;
                }
                if(newMax==0 && commCostList[sstep].rbegin()->a>0)
                {
                    change -= params.L;
                    changing.canShrink = true;
                }
            }
        }

        changing.newCost = cost + change;
        return change;
    }

	// Execute a chosen move, updating the schedule and the data structures
    void executeMove(int node, int newProc, int where, const stepAuxData& changing)
    {
        int oldStep = supstep[node], newStep = oldStep + where;
        int oldProc = proc[node];
        cost = changing.newCost;

        // Work cost change
        workCostList[oldStep].erase(workCostPointer[oldStep][oldProc]);
        workCost[oldStep][oldProc] -= G.workW[node];
        workCostPointer[oldStep][oldProc] = workCostList[oldStep].insert(intPair(workCost[oldStep][oldProc], oldProc)).first;

        workCostList[newStep].erase(workCostPointer[newStep][newProc]);
        workCost[newStep][newProc] += G.workW[node];
        workCostPointer[newStep][newProc] = workCostList[newStep].insert(intPair(workCost[newStep][newProc], newProc)).first;

        // Comm cost change
        for(auto update : changing.sentChange)
            sent[update.first.a][update.first.b] += update.second;
        for(auto update : changing.recChange)
            received[update.first.a][update.first.b] += update.second;

        set<intPair> toUpdate;
        for(auto update : changing.sentChange)
            if(max(sent[update.first.a][update.first.b], received[update.first.a][update.first.b]) != commCost[update.first.a][update.first.b])
                toUpdate.insert(intPair(update.first.a, update.first.b));

        for(auto update : changing.recChange)
            if(max(sent[update.first.a][update.first.b], received[update.first.a][update.first.b]) != commCost[update.first.a][update.first.b])
                toUpdate.insert(intPair(update.first.a, update.first.b));

        for(auto update : toUpdate)
        {
            commCostList[update.a].erase(commCostPointer[update.a][update.b]);
            commCost[update.a][update.b] = max(sent[update.a][update.b], received[update.a][update.b]);
            commCostPointer[update.a][update.b] = commCostList[update.a].insert(intPair(commCost[update.a][update.b], update.b)).first;
        }


        //update successor lists
        for(int pred : G.In[node])
        {
            map<int, int>::iterator itr = succSteps[pred][oldProc].find(oldStep);
            if((--(itr->second))==0)
                    succSteps[pred][oldProc].erase(itr);

            itr = succSteps[pred][newProc].find(newStep);
            if(itr==succSteps[pred][newProc].end())
                    succSteps[pred][newProc].insert({newStep, 1});
                else
                    itr->second += 1;
        }

        //update data
        proc[node] = newProc;
        supstep[node] = newStep;
        supsteplists[oldStep][oldProc].erase(supStepListPointer[node]);
        supsteplists[newStep][newProc].push_back(node);
        supStepListPointer[node] = (--supsteplists[newStep][newProc].end());

        //update movability
        updateNodeMoves(node);
        for(int pred : G.In[node])
            updateNodeMoves(pred);
        for(int succ : G.Out[node])
            updateNodeMoves(succ);

        if(changing.canShrink)
            RemoveNeedlessSupSteps();
    }

	// Single hill climbing step
    bool ImproveBSP(bool greedy = false, bool shrink = true)
    {
        int bestCost = cost;
        stepAuxData bestMoveData;
        intPair bestMove;
        int bestDir;

        for(int dir=0; dir<NumDirections; ++dir)
            for(intPair next : moveOptions[dir])
            {
                stepAuxData moveData;
                int costDiff = moveCostChange(next.a, next.b, dir-1, moveData);

                if(greedy && costDiff<0)
                {
                    executeMove(next.a, next.b, dir-1, moveData);
                    return true;
                }
                else if(cost+costDiff<bestCost)
                {
                    bestCost = cost+costDiff;
                    bestMove = next;
                    bestMoveData = moveData;
                    bestDir = dir-1;
                }
            }


        if(bestCost == cost)
        {
            if(shrink)
            {
                RemoveNeedlessSupSteps();
                return (bestCost != cost);
            }
            else
                return false;
        }

        executeMove(bestMove.a, bestMove.b, bestDir, bestMoveData);


        return true;
    }

	// Main method for hill climbing (with time limit)
    void HillClimb(int TimeLimit = 600, bool greedy = true, bool shrink = true)
    {
        InitHillClimbing();
        steady_clock::time_point startTime = steady_clock::now();

        int counter=0;
        while(ImproveBSP(greedy, shrink))
            if((++counter)==10)
            {
                counter = 0;
                steady_clock::time_point now = steady_clock::now();
                int elapsed = duration_cast<seconds>(now - startTime).count();
                if(elapsed>=TimeLimit)
                {
                    cout<<"Hill Climbing was shut down due to time limit."<<endl;
                    break;
                }
            }

        int computedCost = cost;
        if(computedCost != GetBSPCost())
        {
            cout<<"ERROR: Cost calculation in HillClimbing is incorrect!"<<endl;
        }
    }

	// Hill climbing for limited number of steps
    void HillClimbSteps(int StepsLimit = 10, bool greedy = true, bool shrink = true)
    {
        InitHillClimbing();
        for(int i=0; i<StepsLimit; ++i)
            if(!ImproveBSP(greedy, shrink))
                break;

    }


	// INGREDIENTS OF COMM. SCHEDULE HILL CLIMBING

	// Initialization for comm. schedule hill climbing
    void InitCommSchedule()
    {
        commSchedule.clear();
        commSchedule.resize(G.n, vector<int>(params.p, -1));
        commBounds.clear();
        commBounds.resize(G.n, vector<intPair >(params.p));
        commSchedSendLists.clear();
        commSchedSendLists.resize(supsteplists.size()-1, vector<list<intPair> >(params.p));
        commSchedRecLists.clear();
        commSchedRecLists.resize(supsteplists.size()-1, vector<list<intPair> >(params.p));
        commSchedSendListPointer.clear();
        commSchedSendListPointer.resize(G.n, vector<list<intPair>::iterator>(params.p));
        commSchedRecListPointer.clear();
        commSchedRecListPointer.resize(G.n, vector<list<intPair>::iterator>(params.p));

        for(int i=1; i<supsteplists.size(); ++i)
            for(int j=0; j<params.p; ++j)
                for(int node : supsteplists[i][j])
                    for(int pred : G.In[node])
                        if(proc[pred]!=proc[node] && commSchedule[pred][proc[node]]==-1)
                        {
                            commSchedule[pred][proc[node]] = i-1;
                            commBounds[pred][proc[node]] = intPair(supstep[pred], i-1);

                            commSchedSendLists[i-1][proc[pred]].push_front(intPair(pred, proc[node]));
                            commSchedSendListPointer[pred][proc[node]] = commSchedSendLists[i-1][proc[pred]].begin();
                            commSchedRecLists[i-1][proc[node]].push_front(intPair(pred, proc[node]));
                            commSchedRecListPointer[pred][proc[node]] = commSchedRecLists[i-1][proc[node]].begin();
                        }

    }

	// compute cost change incurred by a potential move
    int CSCostChange(int node, int p, int step)
    {
        int oldStep = commSchedule[node][p];
        int sourceProc = proc[node];
        int change = 0;

        // Change at old place
        auto itr = commCostList[oldStep].rbegin();
        int oldMax = itr->a;
        int maxSource = max(sent[oldStep][sourceProc]-G.commW[node]*params.sendCost[sourceProc][p],received[oldStep][sourceProc]);
        int maxTarget = max(sent[oldStep][p],received[oldStep][p]-G.commW[node]*params.sendCost[sourceProc][p]);
        int maxOther = 0;
        for(; itr!=commCostList[oldStep].rend(); ++itr)
            if(itr->b != sourceProc && itr->b != p)
            {
                maxOther = itr->a;
                break;
            }

        int newMax = max(max(maxSource,maxTarget),maxOther);
        change += newMax - oldMax;

        // Change at new place
        oldMax = commCostList[step].rbegin()->a;
        newMax = max(max(oldMax, sent[step][sourceProc]+G.commW[node]*params.sendCost[sourceProc][p]),received[step][p]+G.commW[node]*params.sendCost[sourceProc][p]);
        change += newMax - oldMax;

        return change;
    }

	// execute a move, updating the comm. schedule and the data structures
    void executeCSMove(int node, int p, int step, int changeCost)
    {
        int oldStep = commSchedule[node][p];
        int sourceProc = proc[node];
        cost += changeCost*params.g;

        // Old step update
        if(sent[oldStep][sourceProc] > received[oldStep][sourceProc])
        {
            commCostList[oldStep].erase(commCostPointer[oldStep][sourceProc]);
            sent[oldStep][sourceProc] -= G.commW[node]*params.sendCost[sourceProc][p];
            commCost[oldStep][sourceProc] = max(sent[oldStep][sourceProc], received[oldStep][sourceProc]);
            commCostPointer[oldStep][sourceProc] = commCostList[oldStep].insert(intPair(commCost[oldStep][sourceProc], sourceProc)).first;
        }
        else
            sent[oldStep][sourceProc] -= G.commW[node]*params.sendCost[sourceProc][p];

        if(received[oldStep][p] > sent[oldStep][p])
        {
            commCostList[oldStep].erase(commCostPointer[oldStep][p]);
            received[oldStep][p] -= G.commW[node]*params.sendCost[sourceProc][p];
            commCost[oldStep][p] = max(sent[oldStep][p], received[oldStep][p]);
            commCostPointer[oldStep][p] = commCostList[oldStep].insert(intPair(commCost[oldStep][p], p)).first;
        }
        else
            received[oldStep][p] -= G.commW[node]*params.sendCost[sourceProc][p];

        //New step update
        sent[step][sourceProc] += G.commW[node]*params.sendCost[sourceProc][p];
        if(sent[step][sourceProc] > received[step][sourceProc])
        {
            commCostList[step].erase(commCostPointer[step][sourceProc]);
            commCost[step][sourceProc] = sent[step][sourceProc];
            commCostPointer[step][sourceProc] = commCostList[step].insert(intPair(commCost[step][sourceProc], sourceProc)).first;
        }

        received[step][p] += G.commW[node]*params.sendCost[sourceProc][p];
        if(received[step][p] > sent[step][p])
        {
            commCostList[step].erase(commCostPointer[step][p]);
            commCost[step][p] = received[step][p];
            commCostPointer[step][p] = commCostList[step].insert(intPair(commCost[step][p], p)).first;
        }

        //CommSched update
        commSchedule[node][p] = step;

        //Comm lists
        commSchedSendLists[oldStep][sourceProc].erase(commSchedSendListPointer[node][p]);
        commSchedSendLists[step][sourceProc].push_front(intPair(node, p));
        commSchedSendListPointer[node][p] = commSchedSendLists[step][sourceProc].begin();

        commSchedRecLists[oldStep][p].erase(commSchedRecListPointer[node][p]);
        commSchedRecLists[step][p].push_front(intPair(node, p));
        commSchedRecListPointer[node][p] = commSchedRecLists[step][p].begin();
    }

	// Single comm. schedule hill climbing step
    bool ImproveCS(bool greedy = false)
    {
        int M = supsteplists.size();
        int bestCost = cost;
        int bestNode, bestProc, bestStep, bestDiff;

        for(int i=0; i<M-1; ++i)
        {
            auto itr = commCostList[i].rbegin();
            int commMax = itr->a;
            if(commMax == 0)
                continue;

            for(; itr != commCostList[i].rend() && itr->a == commMax; ++itr)
            {
                int maxProc = itr->b;
                if(sent[i][maxProc] == commMax)
                    for(intPair entry : commSchedSendLists[i][maxProc])
                    {
                        int node = entry.a;
                        int p = entry.b;
                        for(int step = commBounds[node][p].a; step<commBounds[node][p].b; ++step)
                        {
                            if(step == commSchedule[node][p])
                                continue;

                            int costDiff = CSCostChange(node, p, step);

                            if(greedy && costDiff<0)
                            {
                                executeCSMove(node, p, step, costDiff);

                                return true;
                            }
                            else if(cost+costDiff<bestCost)
                            {
                                bestCost = cost+costDiff;
                                bestNode = node;
                                bestProc = p;
                                bestStep = step;
                                bestDiff = costDiff;
                            }
                        }
                    }

                if(received[i][maxProc] == commMax)
                    for(intPair entry : commSchedRecLists[i][maxProc])
                    {
                        int node = entry.a;
                        int p = entry.b;
                        for(int step = commBounds[node][p].a; step<commBounds[node][p].b; ++step)
                        {
                            if(step == commSchedule[node][p])
                                continue;

                            int costDiff = CSCostChange(node, p, step);

                            if(greedy && costDiff<0)
                            {
                                executeCSMove(node, p, step, costDiff);

                                return true;
                            }
                            else if(cost+costDiff<bestCost)
                            {
                                bestCost = cost+costDiff;
                                bestNode = node;
                                bestProc = p;
                                bestStep = step;
                                bestDiff = costDiff;
                            }
                        }
                    }

            }
        }

        if(bestCost == cost)
            return false;

        executeCSMove(bestNode, bestProc, bestStep, bestDiff);

        return true;
    }

	// Main function for comm. schedule hill climbing
    void CommHillClimb( int TimeLimit = 60, bool greedy = true)
    {
        InitCommSchedule();
        steady_clock::time_point startTime = steady_clock::now();

        int counter=0;
        while(ImproveCS(greedy))
            if((++counter)==1)
            {
                counter = 0;
                steady_clock::time_point now = steady_clock::now();
                int elapsed = duration_cast<seconds>(now - startTime).count();
                if(elapsed>=TimeLimit)
                {
                    cout<<"Comm. Sched. Hill Climbing was shut down due to time limit."<<endl;
                    break;
                }
            }

    }

	// FURTHER AUXILIARY FOR HILL CLIMBING

	// Combine subsequent supersteps whenever there is no communication happening inbetween
    void RemoveNeedlessSupSteps()
    {
        int step = 0;
        if(commSchedule.empty()) // lazy data sending - default comm schedule
        {
            int nextBreak = supsteplists.size();
            for(int i=0; i<supsteplists.size(); ++i)
            {
                if(nextBreak==i)
                {
                    ++step;
                    nextBreak = supsteplists.size();
                }
                for(int j=0; j<params.p; ++j)
                    for(int node : supsteplists[i][j])
                    {
                        supstep[node]=step;
                        for(int succ : G.Out[node])
                            if(proc[node]!=proc[succ] && supstep[succ]<nextBreak)
                                nextBreak = supstep[succ];
                    }
            }
        }
        else //concrete comm schedule
        {
            vector<bool> emptyStep(supsteplists.size(), true);
            for(int i=0; i<G.n; ++i)
                for(int j=0; j<params.p; ++j)
                    if(commSchedule[i][j]>=0)
                        emptyStep[commSchedule[i][j]] = false;

            vector<int> newIdx(supsteplists.size());
            for(int i=0; i<supsteplists.size(); ++i)
            {
                newIdx[i]=step;
                if(!emptyStep[i])
                    ++step;
            }
            for(int i=0; i<G.n; ++i)
            {
                supstep[i]=newIdx[supstep[i]];
                for(int j=0; j<params.p; ++j)
                    if(commSchedule[i][j]>=0)
                        commSchedule[i][j]=newIdx[commSchedule[i][j]];
            }
        }

        // update data structures
        CreateSupStepLists();
        InitHillClimbing();
    }


	// MAIN FUNCTIONS FOR SCHEDULES (R/W, CHECKS)

	// reading scheduel from file
    bool read(string filename, bool isBSP = true, bool NoNUMA = false)
    {
        ifstream infile(filename);
        if(!infile.is_open())
        {
            cout<<"Unable to find/open input schedule file.\n";
            return false;
        }

        G.read(infile);
        readProblemParams(infile, params, NoNUMA);

        int N = G.n;
        proc.clear();
        proc.resize(N);
        time.clear();
        supstep.clear();
        vector<int> entries(N);
        string line;

        //read schedule
        for(int i=0; i<N; ++i)
        {
            if(infile.eof())
            {
                cout<<"Incorrect input file format (file terminated too early).\n";
                return false;
            }

            getline(infile, line);
            while(!infile.eof() && line.at(0)=='%')
                getline(infile, line);

            int node, inProc, inStep;
            sscanf(line.c_str(), "%d %d %d", &node, &inProc, &inStep);

            if(node<0 || inProc<0 || inStep<0 || node>=N || inProc>=params.p)
            {
                cout<<"Incorrect input file format (index out of range for one of the schedule entries)."<<endl;
                return false;
            }

            entries[node] = inStep;
            proc[node] = inProc;
        }
        infile.close();

        if(isBSP)
        {
            supstep = entries;
            CreateSupStepLists();
            if(!CheckBSP())
                return false;
        }
        else
        {
            time = entries;
            ComputeSuperSteps();
            if(!CheckStandard(params.g))
                return false;
        }

        return true;
     }

	// auxiliary for classical schedule validity check
     bool CheckJobOverlap()
     {
        int N=G.n;
        vector<vector <intPair> > jobs(params.p);

        for(int i=0; i<N; ++i)
        {
            intPair next;
            next.a = time[i];
            next.b = next.a + G.workW[i];
        }
        for(int i=0; i<params.p; ++i)
            if(!isDisjoint(jobs[i]))
            {
                cout<<"This is not a valid scheduling (jobs overlap at processor "<<i<<")."<<endl;
                return false;
            }

        return true;
     }

	// check if a classical (non-BSP) schedule is valid
     bool CheckStandard(int delay)
     {
         if(!CheckJobOverlap())
            return false;

         int N = G.n;
         for(int i=0; i<N; ++i)
             for(int j=0; j<G.Out[i].size(); ++j)
             {
                 int node = G.Out[i][j];
                 int diff = (proc[i]==proc[node]) ? 0 : delay*G.commW[i];
                 if(time[i]+G.workW[i]+diff>time[node])
                 {
                     cout<<"This is not a valid scheduling (problems with nodes "<<i<<" and "<<node<<")."<<endl;
                     return false;
                 }
             }

         return true;
     }

	// check if BSP schedule is valid
     bool CheckBSP()
     {
         int N = G.n;
         for(int i=0; i<N; ++i)
         {
             for(int succ : G.Out[i])
             {
                 int diff = (proc[i]==proc[succ]) ? 0 : 1;
                 if(supstep[i]+diff>supstep[succ])
                 {
                     cout<<"This is not a valid scheduling (problems with nodes "<<i<<" and "<<succ<<")."<<endl;
                     return false;
                 }

                 if(!commSchedule.empty() && proc[i]!=proc[succ])
                     if(commSchedule[i][proc[succ]]<supstep[i] || commSchedule[i][proc[succ]]>=supstep[succ])
                     {
                         cout<<"This is not a valid scheduling (problems with nodes "<<i<<" and "<<succ<<")."<<endl;
                         return false;
                     }
             }
         }

         return true;
     }

	// write BSP (problem and) schedule to file
     bool WriteBSP(string filename, bool NoNUMA=false)
     {
        ofstream outfile(filename);
        if(!outfile.is_open())
        {
            cout<<"Unable to write/open output schedule file.\n";
            return false;
        }

        G.write(outfile);

        params.write(outfile, NoNUMA);

        for(int i=0; i<G.n; ++i)
            outfile<<i<<" "<<proc[i]<<" "<<supstep[i]<<endl;

        if(!commSchedule.empty())
        {
            int countLines=0;
            for(int i=0; i<G.n; ++i)
                for(int j=0; j<params.p; ++j)
                    if(commSchedule[i][j]>=0)
                        ++countLines;

            outfile<<countLines<<endl;
            for(int i=0; i<G.n; ++i)
                for(int j=0; j<params.p; ++j)
                    if(commSchedule[i][j]>=0)
                        outfile<<i<<" "<<proc[i]<<" "<<j<<" "<<commSchedule[i][j]<<endl;
        }

        outfile.close();
        return true;
     }



     // COST CALCULATIONS

	// compute classical schedule makespan
     int GetStandardCost()
     {
         int mx = 0;
         for(int i=0; i<G.n; ++i)
            if(time[i]+G.workW[i]> mx)
                mx = time[i]+G.workW[i];

         return mx;
     }

	// compute BSP schedule cost
    int GetBSPCost()
    {
        CheckBSP();

        int nrSupSteps = supsteplists.size();
        int N=G.n;
        int cost=0;

        for(int i=0; i<nrSupSteps; ++i)
        {
            int maxWork = 0;
            for(int j=0; j<params.p; ++j)
            {
                int work=0;
                for(int node : supsteplists[i][j])
                    work += G.workW[node];

                if(work>maxWork)
                    maxWork = work;
            }

            cost+=maxWork;
        }
        if(commSchedule.empty()) // lazy data sending - default comm schedule
        {
            vector<vector <bool> > present(N, vector<bool>(params.p, false));
            for(int i=1; i<nrSupSteps; ++i)
            {
                vector<int> send(params.p, 0), receive(params.p, 0);

                for(int j=0; j<params.p; ++j)
                {
                    for(int target : supsteplists[i][j])
                    {
                        for(int l=0; l<G.In[target].size(); ++l)
                        {
                            int source = G.In[target][l];
                            if(proc[target]!=proc[source] && !present[source][proc[target]])
                            {
                                present[source][proc[target]] = true;
                                send[proc[source]] += G.commW[source] * params.sendCost[proc[source]][proc[target]];
                                receive[proc[target]] += G.commW[source] * params.sendCost[proc[source]][proc[target]];
                            }
                        }
                    }
                }

                int mx=0;
                for(int j=0; j<params.p; ++j)
                {
                    if(send[j]>mx)
                        mx=send[j];
                    if(receive[j]>mx)
                        mx=receive[j];
                }

                int latency = mx > 0 ? params.L : 0;
                cost+=params.g * mx + latency;
            }
        }
        else // we have a specifically defined comm schedule
        {
            vector<vector<int> > send(nrSupSteps-1, vector<int>(params.p, 0)), receive(nrSupSteps-1, vector<int>(params.p, 0));
            for(int i=0; i<N; ++i)
                for(int j=0; j<params.p; ++j)
                    if(commSchedule[i][j]>=0)
                    {
                        send[commSchedule[i][j]][proc[i]] += G.commW[i] * params.sendCost[proc[i]][j];
                        receive[commSchedule[i][j]][j] += G.commW[i] * params.sendCost[proc[i]][j];
                    }

            for(int i=1; i<nrSupSteps; ++i)
            {
                int mx=0;
                for(int j=0; j<params.p; ++j)
                {
                    if(send[i-1][j]>mx)
                        mx=send[i-1][j];
                    if(receive[i-1][j]>mx)
                        mx=receive[i-1][j];
                }

                int latency = mx > 0 ? params.L : 0;
                cost+=params.g * mx + latency;
            }
        }

        return cost;
    }


    // CONVERSIONS

	// create superstep lists (for convenience) for a BSP schedule
    void CreateSupStepLists()
    {
        int N = G.n;
        int nrSupSteps = 0;
        for(int i=0; i<N; ++i)
            if(supstep[i]>=nrSupSteps)
                nrSupSteps = supstep[i]+1;

        supsteplists.clear();
        supsteplists.resize(nrSupSteps, vector<list<int> >(params.p));

        vector<vector<int> > timer(nrSupSteps, vector<int>(params.p, 0));

        vector<int> topOrder = G.GetTopOrder();
        for(int i=0; i<N; ++i)
        {
            int node = topOrder[i];
            supsteplists[supstep[node]][proc[node]].push_back(node);
        }
    }

	// convert a classical schedule (with concrete time points) into a BSP schedule (with supersteps)
    void ComputeSuperSteps()
    {
        int N = G.n;
        supstep.clear();
        supstep.resize(N);

        int stepIdx = 0, totalDone = 0;
        vector<bool> processed(N, false);
        vector<list<int>::iterator> done(params.p), limit(params.p);
        for(int j=0; j<params.p; ++j)
            done[j] = greedyProcLists[j].begin();

        while(totalDone<N)
        {
            int timeLimit = INT_MAX;
            for(int j=0; j<params.p; ++j)
            {
                for(limit[j]=done[j]; limit[j]!=greedyProcLists[j].end(); ++limit[j])
                {
                    int node = *limit[j];
                    bool cut = false;
                    for(int source : G.In[node])
                        if(!processed[source] && proc[source]!=proc[node])
                            cut=true;

                    if(cut)
                        break;
                }
                if(limit[j]!=greedyProcLists[j].end() && time[*limit[j]]<timeLimit)
                    timeLimit = time[*limit[j]];
            }

            for(int j=0; j<params.p; ++j)
                for(; done[j]!=limit[j] && (time[*done[j]]<timeLimit || (time[*done[j]]==timeLimit && G.workW[*done[j]]==0)); ++done[j])
                {
                    processed[*done[j]] = true;
                    supstep[*done[j]]=stepIdx;
                    ++totalDone;
                }

            ++stepIdx;
        }

        CreateSupStepLists();
    }


	// AUXILIARY FUNCTIONS FOR DESIGNING SIMPLE (GREEDY) SCHEDULES

    // Choosing a node to assing for classical greedy methods (cilk, SJF, random)
    void Choose(string mode, const set<int>& readyNodes, const vector<bool>& procFree, int& node, int& p)
    {
        if(mode=="SJF")
        {
            node =-1;
            for(auto& r : readyNodes)
                if(node==-1 || G.workW[r]<G.workW[node])
                    node = r;

            p=0;
            for(; p<params.p; ++p)
                if(procFree[p])
                    break;

        }
        else if(mode=="random")
        {
            int i=0, rnd = randInt(readyNodes.size());
            for(auto& r : readyNodes)
            {
                if(i==rnd)
                {
                    node = r;
                    break;
                }
                ++i;
            }

            int cnt = 0;
            for(int i=0; i<params.p; ++i)
                if(procFree[i])
                    ++cnt;

            rnd = randInt(cnt);
            cnt=0;
            for(int i=0; i<params.p; ++i)
                if(procFree[i])
                {
                    if(cnt==rnd)
                    {
                        p = i;
                        break;
                    }
                    ++cnt;
                }
        }
        else if(mode=="cilk")
        {
            for(int i=0; i<params.p; ++i)
                if(procFree[i] && !procQueue[i].empty())
                {
                    p = i;
                    node = procQueue[i].back();
                    procQueue[i].pop_back();
                    return;
                }

            // Time to steal
            vector<bool> canSteal(params.p, false);
            int cnt = 0;
            for(int i=0; i<params.p; ++i)
            {
                canSteal[i]=!procQueue[i].empty();
                if(canSteal[i])
                    ++cnt;
            }

            for(int i=0; i<params.p; ++i)
                if(procFree[i])
                {
                    p=i;
                    break;
                }

            int rnd = randInt(cnt), q;
            cnt = 0;
            for(int i=0; i<params.p; ++i)
                if(canSteal[i])
                {
                    if(cnt==rnd)
                    {
                        q = i;
                        break;
                    }
                    ++cnt;
                }

            node = procQueue[q].front();
            procQueue[q].pop_front();
        }
    }


    // Choosing a node to assing for BSPg
    void BSPchoose(string mode, const set<int>& allReady, const vector<set<int> > & procReady, const vector<bool>& procFree, int& node, int& p)
    {
		double maxScore = -1;
		for(int i=0; i<params.p; ++i)
			if(procFree[i] && !procReady[i].empty())
			{
				for(auto& r : procReady[i])
				{
					double score = 0;
					for(int j : G.In[r])
						if(procInHyperedge[j][i])
							score += (double)G.commW[j]/(double)G.Out[j].size();


					if(score>maxScore)
					{
						maxScore=score;
						node = r;
						p = i;
					}
				}
				return;
			}

		for(auto& r : allReady)
		{
			for(int i=0; i<params.p; ++i)
			{
				if(!procFree[i])
					continue;

				double score = 0;
				for(int j : G.In[r])
					if(procInHyperedge[j][i])
						score += (double)G.commW[j]/(double)G.Out[j].size();


				if(score>maxScore)
				{
					maxScore=score;
					node = r;
					p = i;
				}
			}
		}
    }

	// auxiliary - check if it is possible to assing a node at all
    bool CanChooseNode(const set<int>& allReady, const vector<set<int> > & procReady, const vector<bool>& procFree)
    {
        for(int i=0; i<params.p; ++i)
            if(procFree[i] && !procReady[i].empty())
                return true;

        if(!allReady.empty())
            for(int i=0; i<params.p; ++i)
                if(procFree[i])
                    return true;

        return false;
    }
};

// INTIALIZER AND BASELINE SHCEDULING METHODS

// Classical greedy scheduler that assigns a new node to processor when previous task is finished
// (covers cilk, SJF and random)
schedule GreedySchedule(const DAG& G, BSPproblem params, string mode)
{
    int N = G.n;
    schedule S;
    S.G = G;
    S.params = params;
    S.proc.clear();
    S.proc.resize(N, -1);
    S.time.clear();
    S.time.resize(N);

    set<int> ready;
    vector<int> predec(N, 0);
    vector<bool> procFree(params.p, true);
    int free = params.p;

    S.greedyProcLists.clear();
    S.greedyProcLists.resize(params.p);

    //aux for comm-aware heuristics
    S.procInHyperedge.clear();
    S.procInHyperedge.resize(N, vector<int>(params.p, false));
    S.procQueue.clear();
    S.procQueue.resize(params.p);

    set<intPair> finishTimes;
    intPair start(0, -1);
    finishTimes.insert(start);

    for(int i=0; i<N; ++i)
    {
        if(G.In[i].empty())
        {
            ready.insert(i);
            if(mode=="cilk")
                S.procQueue[0].push_front(i);
        }
    }

    while(!finishTimes.empty())
    {
        int time = finishTimes.begin()->a;

        // Find new ready jobs
        while(!finishTimes.empty() && finishTimes.begin()->a == time)
        {
            intPair currentPair = *finishTimes.begin();
            finishTimes.erase(finishTimes.begin());
            int node = currentPair.b;
            if(node!=-1)
            {
                for(int j=0; j<G.Out[node].size(); ++j)
                {
                    ++predec[G.Out[node][j]];
                    if(predec[G.Out[node][j]]==G.In[G.Out[node][j]].size())
                    {
                        ready.insert(G.Out[node][j]);
                        if(mode=="cilk")
                            S.procQueue[S.proc[node]].push_back(G.Out[node][j]);
                    }
                }
                procFree[S.proc[node]] = true;
                ++free;
            }

        }


        //Assign new jobs to processors
        while(free>0 && !ready.empty())
        {
            int nextNode, nextProc;
            S.Choose(mode, ready, procFree, nextNode, nextProc);

            ready.erase(nextNode);
            S.proc[nextNode] = nextProc;
            S.time[nextNode] = time;

            finishTimes.insert(intPair(time+G.workW[nextNode], nextNode));
            procFree[nextProc]=false;
            --free;

            // update comm auxiliary structure
            S.greedyProcLists[nextProc].push_back(nextNode);
            S.procInHyperedge[nextNode][nextProc]=true;
            for(int i :G.In[nextNode])
                S.procInHyperedge[i][nextProc]=true;
        }
    }

	// convert the classical schedule to a BSP schedule
    S.ComputeSuperSteps();

    return S;
}

//BSP greedy (BSPg) scheduler
schedule GreedyBSP(const DAG& G, BSPproblem params, string mode)
{
    int N = G.n;
    schedule S;
    S.G = G;
    S.params = params;
    S.proc.clear();
    S.proc.resize(N, -1);
    S.supstep.clear();
    S.supstep.resize(N);

    set<int> ready;

    vector<set<int> > procReady(params.p);
    set<int> allReady;

    vector<int> predec(N, 0);
    vector<bool> procFree(params.p, true);
    int free = params.p;

    S.procInHyperedge.clear();
    S.procInHyperedge.resize(N, vector<int>(params.p, false));

    set<intPair> finishTimes;
    intPair start(0, -1);
    finishTimes.insert(start);

    for(int i=0; i<N; ++i)
        if(G.In[i].empty())
        {
            ready.insert(i);
            allReady.insert(i);
        }

    int supstep = 0;
    bool endSupStep = false;
    while(!ready.empty() || !finishTimes.empty())
    {
        if(finishTimes.empty() && endSupStep)
        {
            for(int i=0; i<params.p; ++i)
                procReady[i].clear();

            allReady = ready;

            ++supstep;
            endSupStep = false;

            intPair start(0, -1);
            finishTimes.insert(start);
        }

        int time = finishTimes.begin()->a;

        // Find new ready jobs
        while(!finishTimes.empty() && finishTimes.begin()->a == time)
        {
            intPair currentPair = *finishTimes.begin();
            finishTimes.erase(finishTimes.begin());
            int node = currentPair.b;
            if(node!=-1)
            {
                for(int j=0; j<G.Out[node].size(); ++j)
                {
                    int succ = G.Out[node][j];
                    ++predec[succ];
                    if(predec[succ]==G.In[succ].size())
                    {
                        ready.insert(succ);

                        bool canAdd = true;
                        for(int i : G.In[succ])
                            if(S.proc[i]!=S.proc[node] && S.supstep[i]==supstep)
                                canAdd= false;

                        if(canAdd)
                            procReady[S.proc[node]].insert(succ);
                    }

                }
                procFree[S.proc[node]] = true;
                ++free;
            }

        }

        if(endSupStep)
            continue;

        //Assign new jobs to processors
        while(true)
        {
            if(!S.CanChooseNode(allReady, procReady, procFree))
                break;

            int nextNode=-1, nextProc;
            S.BSPchoose(mode, allReady, procReady, procFree, nextNode, nextProc);

            if(procReady[nextProc].find(nextNode) != procReady[nextProc].end())
                procReady[nextProc].erase(nextNode);
            else
                allReady.erase(nextNode);

            ready.erase(nextNode);
            S.proc[nextNode] = nextProc;
            S.supstep[nextNode] = supstep;

            finishTimes.insert(intPair(time+G.workW[nextNode], nextNode));
            procFree[nextProc]=false;
            --free;

            // update comm auxiliary structure
            S.procInHyperedge[nextNode][nextProc]=true;
            for(int i :G.In[nextNode])
                S.procInHyperedge[i][nextProc]=true;
        }

        if(allReady.empty() && free>=params.p/2)
            endSupStep = true;

    }

    S.CreateSupStepLists();
    return S;
}

// BASELINE LIST SCHEDULERS: BL-EST AND ETF

// auxiliary: compute bottom level for BL-EST
vector<int> ComputeBottomLevel(const DAG& G, const BSPproblem& params)
{
    vector<int> BL(G.n);
    vector<int> topOrder = G.GetTopOrder();
    for(int i=G.n-1; i>=0; --i)
    {
        int node = topOrder[i];
        int maxval = 0;
        for(int succ : G.Out[node])
            if(BL[succ]>maxval)
                maxval=BL[succ];
        BL[node] = G.workW[node] + maxval * params.avgComm * G.commW[node];
    }
    return BL;
}

// auxiliary: compute EST of node for a given processor
int GetESTforProc(const schedule& S, int node, int proc, int finishTime, const vector<intPair>& predec, vector<int>& send, vector<int>& rec, bool NUMA=false)
{
    int EST = finishTime;
    for(intPair next : predec)
    {
        int t=S.time[next.b]+S.G.workW[next.b];
        if(S.proc[next.b]!=proc)
        {
            t=max(t, send[S.proc[next.b]]);
            t=max(t, rec[proc]);
            t+=S.G.commW[next.b]*(NUMA ? S.params.g*S.params.sendCost[S.proc[next.b]][proc] : S.params.avgComm);
            send[S.proc[next.b]]=t;
            rec[proc]=t;
        }
        EST=max(EST, t);
    }
    return EST;
}

// auxiliary: compute EST of node over all processors
intPair GetBestESTforNode(const schedule& S, int node, const vector<int>& finishTime, const vector<int>& send, const vector<int>& rec, bool NUMA=false)
{
    vector<intPair> predec;
    for(int pred : S.G.In[node])
        predec.push_back(intPair(S.time[pred]+S.G.workW[pred], pred));
    sort(predec.begin(), predec.end());
    vector<int> EST(S.params.p);

    for(int j=0; j<S.params.p; ++j)
    {
        vector<int> newSend=send;
        vector<int> newRec=rec;
        EST[j] = GetESTforProc(S, node, j, finishTime[j], predec, newSend, newRec, NUMA);
    }

    int bestIdx=0;
    for(int j=0; j<S.params.p; ++j)
        if(EST[j]<EST[bestIdx])
            bestIdx=j;

    return intPair(bestIdx, EST[bestIdx]);
}

// run ETF of BL-EST (or their heavily experimental NUMA-aware variants)
schedule RunETF(const DAG& G, const BSPproblem& params, string mode)
{
    schedule S;
    S.G = G;
    S.params = params;
    S.proc.clear();
    S.proc.resize(G.n, -1);
    S.time.clear();
    S.time.resize(G.n);
    vector<int> predecProcessed(G.n, 0);

    S.greedyProcLists.clear();
    S.greedyProcLists.resize(params.p);

    vector<int> finishTimes(params.p, 0), send(params.p, 0), rec(params.p, 0);
    bool NUMA = (mode == "BL-EST-NUMA" || mode == "ETF-NUMA");

    vector<int> BL(G.n, 0);
    if(mode=="BL-EST" or mode=="BL-EST-NUMA")
        BL = ComputeBottomLevel(G, params);

    set<intPair> ready;
    for(int i=0; i<G.n; ++i)
        if(G.In[i].empty())
            ready.insert(intPair(BL[i], i));

    while(!ready.empty())
    {
        int node, proc;
        if(mode=="BL-EST" or mode=="BL-EST-NUMA")
        {
            node = ready.begin()->b;
            ready.erase(ready.begin());
            intPair best = GetBestESTforNode(S, node, finishTimes, send, rec, NUMA);
            proc = best.a;
        }
        if(mode=="ETF" or mode=="ETF-NUMA")
        {
            int bestval = -1;
            for(intPair next : ready)
            {
                int i = next.b;
                vector<intPair> predecList;
                for(int pred : G.In[i])
                    predecList.push_back(intPair(S.time[pred]+S.G.workW[pred], pred));
                sort(predecList.begin(), predecList.end());
                for(int j=0; j<params.p; ++j)
                {
                    vector<int> newSend(params.p);
                    vector<int> newRec(params.p);
                    newSend=send;
                    newRec=rec;
                    int EST = GetESTforProc(S, i, j, finishTimes[j], predecList, newSend, newRec, NUMA);
                    if(bestval==-1 || EST<bestval)
                    {
                        bestval = EST;
                        node=i;
                        proc=j;
                    }
                }
            }
            ready.erase(intPair(0, node));
        }
        S.proc[node] = proc;
        S.greedyProcLists[proc].push_back(node);

        vector<intPair> predecList;
        for(int pred : G.In[node])
            predecList.push_back(intPair(S.time[pred]+S.G.workW[pred], pred));
        sort(predecList.begin(), predecList.end());

        S.time[node] = GetESTforProc(S, node, proc, finishTimes[proc], predecList, send, rec, NUMA);
        finishTimes[proc] = S.time[node]+G.workW[node];

        for(int succ : G.Out[node])
        {
            ++predecProcessed[succ];
            if(predecProcessed[succ]==G.In[succ].size())
                ready.insert(intPair(BL[succ], succ));
        }
    }

    S.ComputeSuperSteps();
    return S;
}

// GENERAL FOR INITIALIZERS

// run the appropriate initializer based on mode
schedule RunGreedyMode(const DAG& G, BSPproblem params, string mode)
{
    if(mode.substr(0,3)=="BSP")
        return GreedyBSP(G, params, mode);
    else if(mode.substr(0,6)=="BL-EST" || mode.substr(0,3)=="ETF")
        return RunETF(G, params, mode);
    else
        return GreedySchedule(G, params, mode);
}

// AUXULIARY FOR RUNNING THE MULTILEVEL ALGORITHM

// given an original DAG G, a schedule on the coarsified G and the contraction steps, project the coarse schedule to the entire G
schedule ComputeUncontractedSchedule(const DAG& G, const schedule& CoarseSchedule, const vector<intPair>& contractionSteps)
{
    vector<int> target = GetFinalImage(G, contractionSteps);

    schedule S;
    S.G = G;
    S.params = CoarseSchedule.params;
    S.proc.clear();
    S.proc.resize(G.n);
    S.supstep.clear();
    S.supstep.resize(G.n);

    for(int i=0; i<G.n; ++i)
    {
        S.proc[i] = CoarseSchedule.proc[target[i]];
        S.supstep[i] = CoarseSchedule.supstep[target[i]];
    }

    if(!CoarseSchedule.commSchedule.empty())
    {
        S.commSchedule.clear();
        S.commSchedule.resize(G.n, vector<int>(S.params.p, -1));
        for(int i=0; i<G.n; ++i)
            for(int j=0; j<S.params.p; ++j)
                if(CoarseSchedule.commSchedule[target[i]][j]>=0)
                {
                    //Check if this comm. step is needed at all after uncoarsening
                    bool needed = false;
                    for(int succ : G.Out[i])
                        if(S.proc[succ] == j)
                        {
                            needed = true;
                            break;
                        }

                    if(needed)
                        S.commSchedule[i][j] = CoarseSchedule.commSchedule[target[i]][j];
                }
    }

    S.CreateSupStepLists();
    return S;
}

// read, compute projected schedule (see above), and write to file
schedule ComputeUncontractedSchedule(const DAG& G, const schedule& CoarseSchedule, string contractFile, string outfilename = "")
{
    vector<intPair> contractionSteps;
    if(!ReadContractionFile(contractFile, contractionSteps))
        return schedule();

    schedule S = ComputeUncontractedSchedule(G, CoarseSchedule, contractionSteps);

    cout<<"Uncontracted cost: "<<S.GetBSPCost()<<endl;

    if(!outfilename.empty())
        S.WriteBSP(outfilename);

    return S;
}

// run refinement: iuncoarsify the DAG in small batches, and apply some steps of hill climbing after each iteration
schedule Refine(const DAG& G, const schedule& CoarseSchedule, string contractFile, int batchLength = 20, int HCsteps = 20)
{
    vector<intPair> contractionSteps;
    if(!ReadContractionFile(contractFile, contractionSteps))
        return schedule();

    schedule S = CoarseSchedule;
    S.commSchedule.clear();

    schedule uncontractedS = ComputeUncontractedSchedule(G, S, contractionSteps);

    for(int nextN = CoarseSchedule.G.n+batchLength; nextN<=G.n; nextN += batchLength)
    {
        if(nextN+batchLength>G.n)
            nextN = G.n;

        vector<intPair> contractionPrefix = contractionSteps;
        contractionPrefix.resize(G.n-nextN);

        DAG nextG = G.Contract(contractionPrefix);
        S.G = nextG;
        S.params = CoarseSchedule.params;
        S.proc.clear();
        S.proc.resize(nextN);
        S.supstep.clear();
        S.supstep.resize(nextN);

        //Project full schedule to current graph
        vector<int> target = GetFinalImage(G, contractionPrefix);
        for(int i=0; i<G.n; ++i)
        {
            S.proc[target[i]] = uncontractedS.proc[i];
            S.supstep[target[i]] = uncontractedS.supstep[i];
        }
        S.CreateSupStepLists();

        S.HillClimbSteps(HCsteps);
        uncontractedS = ComputeUncontractedSchedule(G, S, contractionPrefix);
    }

    cout<<"Refined cost: "<<S.GetBSPCost()<<endl;
    return S;
}

// read problem, coarsify DAG and write to file
void RunContract(string filename, ProgramParams options)
{
    DAG G;
    BSPproblem params(4,1,5);
    if(!readProblem(filename, G, params, options.NoNUMA))
        return;

    CoarsifyAndWrite(G, params, max(min(G.n, options.newN), 10), "coarse_"+filename, "contraction_"+filename);
}

// read problem, refine schedule, run comm. schedule hill climbing and write to file
void RunRefine(string filename, string schedulefile, string contractfile, ProgramParams options)
{
    DAG G;
    BSPproblem params(4,1,5);
    if(!readProblem(filename, G, params, options.NoNUMA))
        return;

    schedule S;
    if(!S.read(schedulefile, true, options.NoNUMA))
        return;

    schedule S2 = Refine(G, S, contractfile, 5, options.HCsteps);

    S2.CommHillClimb(options.TimeLimit);
    if(!options.NoOutfile)
        S2.WriteBSP(editFilename(filename, "REF_"), options.NoNUMA);
}


// MAIN FUNCTIONS TO RUN THE SCHEDULERS

// run program, applying only hill climbing
void RunOnlyHC(string filename, ProgramParams options)
{
    schedule S;
    if(!S.read(filename, true, options.NoNUMA))
        return;

    int baseCost = S.GetBSPCost();
    S.HillClimb(options.TimeLimit*0.9);
    int midCost = S.GetBSPCost();
    if(!options.NoOutfile)
        S.WriteBSP(editFilename(filename, "HC"), options.NoNUMA);

    S.CommHillClimb(options.TimeLimit*0.1);
    if(!options.NoOutfile)
        S.WriteBSP(editFilename(filename, "HCCS"), options.NoNUMA);

    if(options.display)
        cout<<baseCost<< " -> "<<midCost<<" -> "<<S.GetBSPCost()<<endl;
}

// run program, applying the specified initializers and (potentially) hill climbing
void RunAlgorithms(string filename, ProgramParams options)
{
    DAG G;
    BSPproblem params(4,1,5);
    if(!readProblem(filename, G, params, options.NoNUMA))
        return;

    for(string mode : possibleModes)
    {
        if(!options.runMode[mode])
            continue;

        schedule S = RunGreedyMode(G, params, mode);
        if(!options.NoOutfile)
            S.WriteBSP(mode+"_"+filename, options.NoNUMA);

        if(options.NoHC)
        {
            if(options.display)
                cout<<mode<<": "<<S.GetBSPCost()<<endl;
            continue;
        }

        int baseCost = S.GetBSPCost();
        int base_steps = S.supsteplists.size();
        S.HillClimb(options.TimeLimit*0.9);
        int midCost = S.GetBSPCost();
        int HC_steps = S.supsteplists.size();
        if(!options.NoOutfile)
            S.WriteBSP(mode+"HC_"+filename, options.NoNUMA);

        S.CommHillClimb(options.TimeLimit*0.1);
        if(!options.NoOutfile)
            S.WriteBSP(mode+"HCCS_"+filename, options.NoNUMA);

        if(options.display)
            cout<<mode<<": "<<baseCost<< " -> "<<midCost<<" -> "<<S.GetBSPCost()<<"   ["<<base_steps<<" -> "<<HC_steps<<"]"<<endl;

    }
}

// invoked upon program call
int main(int argc, char* argv[])
{
    string infile;
    ProgramParams options;

    string schedulefile, contractfile; //for refining
    set<string> needsValue({"-input", "-timeLimit", "-shrinkTo", "-scheduleFile", "-contractFile", "-batchLength", "-HCsteps"});

    // PROCESS COMMAND LINE ARGUMENTS
    for (int i = 1; i < argc; ++i)
    {
        // Check parameters that require an argument afterwards
        if (needsValue.count(argv[i])==1 && i + 1 >= argc)
        {
            cerr << "Parameter error: no parameter value after the \""<<string(argv[i])<<"\" option." << endl;
            return 1;
        }

        else if (string(argv[i]) == "-input")
            infile = argv[++i];

        else if (string(argv[i]) == "-timeLimit")
            options.TimeLimit = stoi(argv[++i]);

        else if (string(argv[i]) == "-noHC")
            options.NoHC = true;

        else if (string(argv[i]) == "-noOutfile")
            options.NoOutfile = true;

        else if (string(argv[i]) == "-noNUMA")
            options.NoNUMA = true;

        else if (string(argv[i]) == "-display")
            options.display = true;

        else if (string(argv[i]) == "-onlyHC")
            options.OnlyHC = true;

        else if (string(argv[i]) == "-contract")
            options.ContractMode = true;

        else if (string(argv[i]) == "-refine")
            options.RefineMode = true;

        else if (string(argv[i]) == "-shrinkTo")
            options.newN = stoi(argv[++i]);

        else if (string(argv[i]) == "-scheduleFile")
            schedulefile = argv[++i];

        else if (string(argv[i]) == "-contractFile")
            contractfile = argv[++i];

        else if (string(argv[i]) == "-batchLength")
            options.batchLength = stoi(argv[++i]);

        else if (string(argv[i]) == "-HCsteps")
            options.HCsteps = stoi(argv[++i]);

        else
        {
            string s = string(argv[i]);
            if(s.substr(0,1)=="-")
            {
                string mode = s.substr(1, s.size()-1);
                if(find(possibleModes.begin(), possibleModes.end(), mode)!=possibleModes.end())
                {
                    options.runMode[mode] = true;
                    continue;
                }
            }
            cerr << "Parameter error: unknown parameter/option "<< string(argv[i]) << endl;
            return 1;
        }
    }

    if(infile.empty())
    {
        cerr << "Input file name not specified."<< endl;
            return 1;
    }

    bool noMode = true;
    for(auto entry : options.runMode)
        if(entry.second)
            noMode=false;

    if(!options.OnlyHC && !options.ContractMode && !options.RefineMode && noMode)
    {
        cerr << "Parameter error: no scheduling algorithms were specified."<< endl;
        return 1;
    }

    if(options.OnlyHC && options.NoHC)
    {
        cerr << "Parameter error: cannot have noHC and onlyHC together."<< endl;
        return 1;
    }

    if((!noMode || options.OnlyHC)  && (options.ContractMode || options.RefineMode))
    {
        cerr << "Parameter error: no scheduling algorithms can be specified when contracting/refining."<< endl;
        return 1;
    }

    if(options.newN != 0  && !options.ContractMode)
    {
        cerr << "Parameter error: cannot use shrinkTo parameter in non-contract mode."<< endl;
        return 1;
    }

    if(options.newN <= 0 && options.ContractMode)
    {
        cerr << "Parameter error: shrinkTo parameter not specified or invalid."<< endl;
        return 1;
    }
    if(!options.RefineMode && (options.batchLength != 10 || options.HCsteps != 50))
    {
        cerr << "Parameter error: can only use batchLength and HCsteps parameters in refine mode."<< endl;
        return 1;
    }
    if(options.RefineMode && (schedulefile.empty() || contractfile.empty()))
    {
        cerr << "Parameter error: scheduleFile and contractFile parameters are required in refine mode."<< endl;
        return 1;
    }

    if(options.OnlyHC && !noMode)
        cout << "Note: algorithm parameters ignored due to onlyHC mode."<< endl;



    if(options.OnlyHC)
        RunOnlyHC(infile, options);
    else if(options.ContractMode)
        RunContract(infile, options);
    else if(options.RefineMode)
        RunRefine(infile, schedulefile, contractfile, options);
    else
        RunAlgorithms(infile, options);

    return 0;
}
