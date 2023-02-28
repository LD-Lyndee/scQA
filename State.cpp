
#include "State.h"
int NumGene;
int NumCell;
int NumState;
vector<vector<float> >InData; 
vector<vector<int> >RecordState;

void Load(char* Fname, int a) {
    vector<vector<float> >data;
    ifstream f;
    f.open(Fname);
    if (f.is_open() == 0) {
        printf("Fail to open data\n");
        exit(1);
    }
    string str;
    while (getline(f, str))
    {
        istringstream input(str);
        vector<float> tmp;
        float a;
        while (input >> a)
            tmp.push_back(a);
        data.push_back(tmp);
    }
    f.close();
    NumCell = data.size();
    NumGene = data[0].size();
    for (int i = 0; i < NumGene; i++) {
        InData.push_back(vector<float>());
        for (int j = 0; j < NumCell; j++) {
            InData[i].push_back(data[j][i]);
        }
    }
    //initial number of states, default is 6
    NumState = a;
}

void Normalize() {
    for (int i = 0; i < NumCell; i++) {
        float max = 0;
        float min = 100;
        for (int j = 0; j < NumGene; j++) {
            if (InData[j][i] > max)
                max = InData[j][i];
            if (InData[j][i] < min)
                min = InData[j][i];
        }
        for (int j = 0; j < NumGene; j++) {
            InData[j][i] = (InData[j][i] - min) / (max - min);
        }
    }
}

vector<int> DefineState() {
    vector<vector<int> >record(NumGene, vector<int>(NumState, 0));
    float gap = (float)1 / (NumState - 1);
    for (int i = 0; i < NumGene; i++) {
        RecordState.push_back(vector<int>());
        for (int j = 0; j < NumCell; j++) {
            if (InData[i][j] == 0) {
                record[i][0]++;
                RecordState[i].push_back(0);
            }
            for (int p = 1; p <= NumState - 1; p++) {
                if (InData[i][j] > gap * (p - 1) && InData[i][j] <= gap * p) {
                    record[i][p]++;
                    RecordState[i].push_back(p);
                }
            }
        }
    }   
    //delete state if the number of cells is below 0.02*NumCell
    vector<int>FlagState(NumGene, 0);
    for (int i = 0; i < NumGene; i++) {
        for (int j = 0; j < NumState; j++) {
            if (record[i][j] < 0.02 * NumCell) {
                record[i][j] = 0;
            }
            else
                FlagState[i]++;
        }
    }
    return FlagState;
}

vector<int>State(char* Fname, int a) {
    Load(Fname, a);
    Normalize();
    return DefineState();
}


