#pragma once
#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include <sstream>
#include<map>
#include<queue>
#include<math.h>
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS 

using namespace std;
extern int NumCell;
extern int NumGene;
extern vector<vector<float> >InData;
extern int NumState;

class BG
{
private:
    //variable
    vector<vector<int> >BData;
    vector<int>TO;
    vector<vector<float> >Transfer;//gene transfer rate
    vector<vector<float> >TransferC;//cluster transfer rate
    int NumBGene;//number of genes
    int NumBCell;//number of cells
    int GeneClusters;//number of clusters of genes
    vector<float>Per;//%
    vector<vector<float> >Sign;//significant score
    vector<vector<int> > Group;//partition genes*clusters
    vector<vector<int> >Template;//template vector clusters*cells
    vector<int>ClusterIndex;//cluster index of genes
    vector<int>CountGene;//number of genes in the cluster
    vector<int>X;
    vector<bool>MarkedClusters;
    float Rho;//vigil paremeter
    int RESET;//reset parameter
    vector<vector<int> >Record;//number of ones
    vector<vector<int> >BLM;
    int MaxIterate;
    int MinGenes;

    //function
    void Load(vector<vector<float> >&, vector<int>&);
    void Sort();
    vector<vector<float> > DisScore(char&);
    vector<vector<float> > PreMat(vector<vector<float> >&);
    void Initial();
    void RunGene(int);
    int    RunCom(int);
    void Train(int, int);
    int    RunComS(int);
    void TrainS(int, int);
    void WeightedRho(int);
    void DFS(vector<vector<int> >&, vector<int>&, int, int);
    void visit(vector<int>&, int, int);
    void SC();//update sign
    void Second();
    void Prune();
public:
    BG();                           
    vector<vector<int> > QualitativeCluster(vector<vector<float> >&, vector<int>&, char&);
};
