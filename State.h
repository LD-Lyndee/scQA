#pragma once
#ifndef STATE_H_ 
#define STATE_H_
#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include <sstream>

#pragma warning(disable:4996)
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS //strcpy、strcat 等函数不安全，忽略警告

using namespace std;

extern int NumGene;
extern int NumCell;
extern int NumState;
extern vector<vector<float> >InData;
void Load(char*, int);
vector<int>State(char*, int);
#endif


