#pragma once

#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <list>
#include<queue>

#pragma warning(disable:4996)
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS //strcpy、strcat 等函数不安全，忽略警告


using namespace std;

class CC {
private:
	int NumCell;
	int NumB;
	int NumQ;
	vector<vector<int> >BLM;
	vector<vector<float> >QLM;
	vector<vector<float> >G2;
	vector<vector<float> >g2;
	vector<vector<int> >N1;
	vector<vector<int> >N2;
	void Load(vector<vector<int> >&, vector<vector<float> >&);
	vector<vector<float> >SimScore(char&, int);
	void SimGraph();
	void PA();
 
public:
	CC();
	void Run(vector<vector<int> >&, vector<vector<float> >&);

};

