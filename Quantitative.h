#pragma once

#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include <sstream>
#include <map>
#include <math.h>
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS //strcpy、strcat 等函数不安全，忽略警告


using namespace std;


class QG
{
private:
    //variable
    vector<vector<float> >QData;
    int NumQGene;//number of genes
    int NumQCell;//number of cells
    int GeneClusters;//number of clusters of genes
    vector<vector<int> > Group;//partition genes*clusters
    vector<vector<float> >Template;//template vector clusters*cells
    vector<vector<float> >Template2;//template vector clusters*cells
    vector<vector<float> >Template3;//template vector clusters*cells
    vector<int>ClusterIndex;//cluster index of genes
    vector<int>CountGene;//number of genes in the cluster
    vector<float>X;//
    vector<bool>MarkedClusters;
    float Rho;//vigil paremeter
    int RESET;//reset parameter
    vector<vector<vector<float> > >Ssort;
    int MaxIterate;
    int MinGenes;

    //function
    vector<vector<float> > PreMat(vector<vector<float> >&);
    void Load(vector<vector<float> >&, vector<int>&, vector<vector<int> >&);
    void Initial();
    vector<vector<float> > SimilarityScore(char&);//calcualte similarity score
    void RunGene(int);
    int    RunCom(int);//run compare
    void Train(int, int);
    void WeightedRho(int);
    void Generate(int);
    void visit(vector<int>&, int, int);
    void DFS(vector<vector<int> >&, vector<int>&, int, int);
    void Prune();
public:
    QG();                          
    vector<vector<float> > QuantitativeCluster(vector<vector<float> >&, vector<int>&, char&, vector<vector<int> >&);
};

