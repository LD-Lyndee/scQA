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
#define _CRT_SECURE_NO_WARNINGS //strcpy��strcat �Ⱥ�������ȫ�����Ծ���

using namespace std;

extern int NumGene;
extern int NumCell;
extern int NumState;
extern vector<vector<float> >InData;
void Load(char*, int);
vector<int>State(char*, int);
#endif


