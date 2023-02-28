#include"State.h"
#include"Qualitative.h"
#include"Quantitative.h"
#include"Cluster.h"

extern int NumGene;
extern int NumCell;
extern int NumState;
extern vector<vector<float> >InData;
extern vector<vector<int> >RecordState;
int main(int argc, char* argv[]) {
    if (argc > 2) {
        int a = 0;
        if (argv[2][0] != 'o') {
            a = strtod(argv[2], NULL);
        }
        else
            a = 6;
        //gene state
        vector<int>flag = State(argv[1], a);
        vector<vector<int> >Gene(NumState + 1);
        for (int i = 0; i < NumGene; i++) {
            Gene[flag[i]].push_back(i);
        }
        vector<int>temp;
        for (int i = 0; i < NumGene; i++) {
            if (flag[i] >= 2) {
                temp.push_back(i);
            }
        }
        //qualitative cluster of genes
        BG B;
        char m1 = 'M';//MCRFPJ 
        vector<vector<int> >BLM = B.QualitativeCluster(InData, temp, m1);
        //quantitative cluster of genes
        char m = 'S';//MSPKL
        vector<vector<float> >Q(NumCell, vector<float>());
        for (int i = 3; i < Gene.size(); i++) {
            QG A;
            if (Gene[i].size() == 0) {
                continue;
            }
            vector<vector<float> >tt = A.QuantitativeCluster(InData, Gene[i], m, RecordState);
            for (int j = 0; j < NumCell; j++) {
                for (int p = 0; p < tt[j].size(); p++) {
                    Q[j].push_back(tt[j][p]);
                }
            }        
        }
        //cluster cells
        CC C;
        C.Run(BLM, Q);
    }
    else {
        printf("MISSING_PARAMETER\n");
        exit(0);
    }
}
