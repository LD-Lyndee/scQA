
//qualitative cluster of genes
//METHOD DEFINITIONS

#include"Qualitative.h"

BG::BG() {
	RESET = 1;
	GeneClusters = 0;
	Rho = 0.7;
	MaxIterate = 1;
	MinGenes = 3;//least number of genes in a cluster 3
}

//preprocessing matrix, normalize to 0-1
vector<vector<float> >BG::PreMat(vector<vector<float> >& data) {
	vector<vector<float> >temp;
	for (int i = 0; i < data.size(); i++) {
		temp.push_back(vector<float>());
		float max = 0;
		for (int j = 0; j < data[i].size(); j++) {
			if (data[i][j] > max) {
				max = data[i][j];
			}
		}
		for (int j = 0; j < data[i].size(); j++) {
			temp[i].push_back(data[i][j] / max);
		}
	}
	return temp;
}

//load two-state genes
void BG::Load(vector<vector<float> >& data, vector<int>& flag) {
	//set filter
	float f = 0; 
	for (int i = 0; i < flag.size(); i++) {
		BData.push_back(vector<int>());
		Transfer.push_back(vector<float>());
		for (int j = 0; j < data[0].size(); j++) {
			Transfer[i].push_back(data[flag[i]][j]);
			if (data[flag[i]][j] <= f) {
				BData[i].push_back(0);
			}
			else {
				BData[i].push_back(1);

			}
		}
	}
	NumBGene = flag.size();
	NumBCell = data[0].size();
}

void BG::Initial() {
	for (int i = 0; i < NumBGene; i++) {
		ClusterIndex.push_back(-1);
	}
	for (int i = 0; i < NumBCell; i++) {
		Per.push_back(0);
	}
	for (int i = 0; i < NumBGene; i++) {
		Group.push_back(vector<int>());
	}
}

void BG::WeightedRho(int j) {
	Rho = Rho - j * 0.05;
	Rho = max(float(0.7), Rho);
}

void BG::Sort() {
	//sort genes by number of non-zero
	vector<int>ss;
	vector<pair<int, int> >ps;
	for (int i = 0; i < NumBGene; i++) {
		ss.push_back(0);
		for (int j = 0; j < NumBCell; j++) {
			if (BData[i][j] == 1) {
				ss[i]++;
				Per[j]++;
			}
		}
		ps.push_back(make_pair(ss[i], i));
	}
	for (int i = 0; i < NumBCell; i++) {
		Per[i] = Per[i] / NumBGene;
	}
	//sort
	sort(ps.begin(), ps.end());
	vector<vector<int> >temp;
	for (int i = 0; i < NumBGene; i++) {
		temp.push_back(vector<int>());
		TO.push_back(ps[i].second);
		for (int j = 0; j < NumBCell; j++) {
			temp[i].push_back(BData[ps[i].second][j]);
		}
	}
	BData = temp;
}

vector<vector<float> > BG::DisScore(char& m) {
	//distance matrix
	//only calculate similarity between neighbors
	int cn = 0.05 * NumBGene;//0.1,0.2
	vector<float>scores(NumBGene, 0);
	vector<vector<float> >score(NumBGene, vector<float>(NumBGene, 0));
	switch (m) {
	case'M':
		//score with Hamming Distance
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < i + cn; j++) {
				if (j >= NumBGene) {
					break;
				}
				for (int p = 0; p < NumBCell; p++) {
					score[i][j] = score[i][j] + (1 - fabs(BData[i][p] - BData[j][p]));
				}
				score[j][i] = score[i][j];
			}
		}
		break;
	case'C':
		//score with chi-square
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < i + cn; j++) {
				if (j >= NumBGene) {
					break;
				}
				float a = 0, b = 0, c = 0, d = 0;
				for (int p = 0; p < NumBCell; p++) {
					if (BData[i][p] == 1 && BData[j][p] == 1) {
						a++;
					}
					else if (BData[i][p] == 0 && BData[j][p] == 1) {
						b++;
					}
					else if (BData[i][p] == 1 && BData[j][p] == 0) {
						c++;
					}
					else
						d++;
				}
				score[i][j] = ((a * d - b * c) / (a + b) * (a + c)) * ((a * d - b * c) / (b + d) * (c + d));
				score[j][i] = score[i][j];
			}
		}
		break;
	case'R':
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < i + cn; j++) {
				if (j >= NumBGene) {
					break;
				}
				int a = 0, b = 0;
				for (int p = 0; p < NumBCell; p++) {
					if (BData[i][p] == BData[j][p]) {
						a++;
					}
					else
						b++;
				}
				score[i][j] = (float)a / (a + b);
				score[j][i] = score[i][j];
			}
		}
		break;
	case'F':
		//score with F1-score
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < i + cn; j++) {
				if (j >= NumBGene) {
					break;
				}
				int a = 0, b = 0, c = 0, d = 0;
				for (int p = 0; p < NumBCell; p++) {
					if (BData[i][p] == 1 && BData[j][p] == 1) {
						a++;
					}
					else if (BData[i][p] == 0 && BData[j][p] == 1) {
						b++;
					}
					else if (BData[i][p] == 1 && BData[j][p] == 0) {
						c++;
					}
					else
						d++;
				}
				score[i][j] = (float)(2 * a) / (float)(2 * a + b + c);
				score[j][i] = score[i][j];
			}
		}
		break;
	case'P':
		//score with correlation 
		for (int i = 0; i < NumBGene; i++) {
			for (int j = 0; j < NumBCell; j++) {
				scores[i] = scores[i] + BData[i][j];
			}
			scores[i] = scores[i] / NumBCell;
		}
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < i + cn; j++) {
				if (j >= NumBGene) {
					break;
				}
				for (int p = 0; p < NumBCell; p++) {
					if (BData[i][p] == 1 && BData[j][p] == 1) {
						score[i][j] = score[i][j] + 1;
					}
				}
				score[i][j] = score[i][j] - NumBCell * scores[i] * scores[j];
				score[i][j] = score[i][j] / (sqrt((scores[i] - scores[i] * scores[i]) * (scores[j] - scores[j] * scores[j])));
				score[j][i] = score[i][j];
			}
		}
		break;
	case'J':
		//score with  cosine
		for (int i = 0; i < NumBGene; i++) {
			for (int j = 0; j < NumBCell; j++) {
				scores[i] = scores[i] + BData[i][j];
			}
			scores[i] = sqrt(scores[i]);
		}
		for (int i = 0; i < NumBGene - 1; i++) {
			for (int j = i + 1; j < NumBGene; j++) {
				for (int p = 0; p < NumBCell; p++) {
					if (BData[i][p] == 1 && BData[j][p] == 1) {
						score[i][j] = score[i][j] + 1;
					}
				}
				score[i][j] = score[i][j] / (scores[i] * scores[j]);
				score[j][i] = score[i][j];
			}
		}
		break;
	}
	return score;
}

void BG::RunGene(int tp) {
	X = BData[tp];
	int j = RunCom(tp);
	Train(tp, j);
	X.clear();
	RESET = 1;
}

int BG::RunCom(int tp) {
	//add transfer rate
	if (GeneClusters == 0) {
		Template.push_back(vector<int>());
		TransferC.push_back(vector<float>());
		for (int i = 0; i < NumBCell; i++) {
			Template[0].push_back(X[i]);
			TransferC[0].push_back(Transfer[tp][i]);
		}
		for (int i = 0; i < NumBGene; i++) {
			Group[i].push_back(0);
		}
		Group[tp][0] = 1;
		ClusterIndex[tp] = GeneClusters;
		Record.push_back(vector<int>());
		for (int i = 0; i < NumBCell; i++) {
			Record[GeneClusters].push_back(BData[tp][i]);
		}
		MarkedClusters.push_back(true);
		CountGene.push_back(1);
		GeneClusters++;
	}
	else {
		vector<vector<float> >Choice(2);
		int index = 0;
		for (int i = 0; i < GeneClusters; i++) {
			if (MarkedClusters[i]) {
				Choice[1].push_back(i);
				Choice[0].push_back(0);
				float temp = 0;
				for (int j = 0; j < NumBCell; j++) {
					if ((X[j] == 1 && Template[i][j] == 1) || (X[j] == 0 && Template[i][j] == 0)) {
						Choice[0][index] = Choice[0][index] + 1;
					}
					else if (X[j] == 1 && Template[i][j] == 0) {
					}
					else if (X[j] == 0 && Template[i][j] == 1) {
						if (TransferC[i][j] < 0.5) {
							Choice[0][index] = Choice[0][index] + TransferC[i][j];
						}
					}
				}
				Choice[0][index] = Choice[0][index] / NumBCell;
				index++;
			}
		}
		while (RESET) {
			float max = Rho;
			int maxj = -1;
			int maxi = -1;
			for (int i = 0; i < Choice[1].size(); i++) {
				if (Choice[0][i] > max) {
					max = Choice[0][i];
					maxi = i;
					maxj = Choice[1][i];
				}
			}
			if (maxj != -1) {
				RESET = 0; //resonate
				return maxj;
			}
			else {
				RESET = -1;//new cluster
				return maxj;
			}
		}
	}
}

void BG::Train(int tp, int j) {
	//update law
	if (RESET == 0) {
		//resonate
		Group[tp][j] = 1;
		ClusterIndex[tp] = j;
		for (int i = 0; i < NumBCell; i++) {
			Record[j][i] = Record[j][i] + BData[tp][i];
		}
		//update template
		CountGene[j]++;
		for (int i = 0; i < NumBCell; i++) {
			if ((CountGene[j] - Record[j][i]) > Record[j][i]) {
				Template[j][i] = 0;
			}
			else {
				Template[j][i] = 1;
			}
			if (Transfer[tp][i] == 1) {
				TransferC[j][i] = (TransferC[j][i] * (Record[j][i] - 1) + Transfer[tp][i]) / Record[j][i];
			}
		}
	}
	else {
		CountGene.push_back(1);
		Template.push_back(vector<int>());
		TransferC.push_back(vector<float>());
		for (int i = 0; i < NumBCell; i++) {
			Template[GeneClusters].push_back(BData[tp][i]);
			TransferC[GeneClusters].push_back(Transfer[tp][i]);
		}
		for (int i = 0; i < NumBGene; i++) {
			Group[i].push_back(0);
		}
		Group[tp][GeneClusters] = 1;
		ClusterIndex[tp] = GeneClusters;
		Record.push_back(vector<int>());
		for (int i = 0; i < NumBCell; i++) {
			Record[GeneClusters].push_back(BData[tp][i]);
		}
		MarkedClusters.push_back(true);
		GeneClusters++;
	}
}

int BG::RunComS(int tp) {
	//add transfer rate
	vector<vector<float> >Choice(2);
	int index = 0;
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			Choice[1].push_back(i);
			Choice[0].push_back(0);
			float temp = 0;
			for (int j = 0; j < NumBCell; j++) {
				if ((X[j] == 1 && Template[i][j] == 1) || (X[j] == 0 && Template[i][j] == 0)) {
					Choice[0][index] = Choice[0][index] + Sign[i][j];
					temp = temp + Sign[i][j];
				}
				else if (X[j] == 1 && Template[i][j] == 0) {
					temp = temp + Sign[i][j];
				}
				else if (X[j] == 0 && Template[i][j] == 1) {
					temp = temp + Sign[i][j];
					if (TransferC[i][j] < 0.5) {
						Choice[0][index] = Choice[0][index] + TransferC[i][j] * Sign[i][j];
					}
				}
			}
			Choice[0][index] = Choice[0][index] / temp;
			index++;
		}
	}
	while (RESET) {
		multimap<float, int>temp;
		for (int i = 0; i < Choice[1].size(); i++) {
			temp.insert(make_pair(Choice[0][i], Choice[1][i]));
		}
		multimap<float, int>::reverse_iterator iter = temp.rbegin();
		multimap<float, int>::reverse_iterator iter2 = temp.rbegin()++;
		if (iter->first < Rho) {
			RESET = -1;
			return RESET;
		}
		else if (iter2->first > Rho) {
			RESET = 1;
			//resonate
			Group[tp][iter->second] = 1;
			Group[tp][iter2->second] = 1;
			ClusterIndex[tp] = iter->second;
			for (int i = 0; i < NumBCell; i++) {
				Record[iter->second][i] = Record[iter->second][i] + BData[tp][i];
			}
			//update template
			CountGene[iter->second]++;
			for (int i = 0; i < NumBCell; i++) {
				Record[iter2->second][i] = Record[iter2->second][i] + BData[tp][i];
			}
			//update template
			CountGene[iter2->second]++;
			for (int i = 0; i < NumBCell; i++) {
				if ((CountGene[iter->second] - Record[iter->second][i]) > Record[iter->second][i]) {
					Template[iter->second][i] = 0;
				}
				else {
					Template[iter->second][i] = 1;
				}
				if (Transfer[tp][i] == 1) {
					TransferC[iter->second][i] = (TransferC[iter->second][i] * (Record[iter->second][i] - 1) + Transfer[tp][i]) / Record[iter->second][i];
				}
			}
			for (int i = 0; i < NumBCell; i++) {
				if ((CountGene[iter2->second] - Record[iter2->second][i]) > Record[iter2->second][i]) {
					Template[iter2->second][i] = 0;
				}
				else {
					Template[iter2->second][i] = 1;
				}
				if (Transfer[tp][i] == 1) {
					TransferC[iter2->second][i] = (TransferC[iter2->second][i] * (Record[iter2->second][i] - 1) + Transfer[tp][i]) / Record[iter2->second][i];
				}
			}
			return RESET;
		}
		else {
			RESET = 0;
			return iter->second;
		}
	}
}

void BG::TrainS(int tp, int j) {
	//update law
	//number of 0 larger than twice number 1
	if (RESET == 0) {
		//resonate
		Group[tp][j] = 1;
		ClusterIndex[tp] = j;
		for (int i = 0; i < NumBCell; i++) {
			Record[j][i] = Record[j][i] + BData[tp][i];
		}
		//update template
		CountGene[j]++;
		for (int i = 0; i < NumBCell; i++) {
			if ((CountGene[j] - Record[j][i]) > Record[j][i]) {
				Template[j][i] = 0;
			}
			else {
				Template[j][i] = 1;
			}
			if (Transfer[tp][i] == 1) {
				TransferC[j][i] = (TransferC[j][i] * (Record[j][i] - 1) + Transfer[tp][i]) / Record[j][i];
			}
		}
	}
}

void BG::Prune() {
	//Prune clusters
	for (int i = 0; i < GeneClusters; i++) {
		if (CountGene[i] < MinGenes) {
			MarkedClusters[i] = false;
			for (int j = 0; j < NumBGene; j++) {
				if (ClusterIndex[j] == i) {
					ClusterIndex[j] = -1;
					Group[j][i] = 0;
				}
			}
		}
	}
}

vector<vector<int> > BG::QualitativeCluster(vector<vector<float> >& data, vector<int>& flag, char& m) {
	vector<vector<float> >data2 = PreMat(data);
	Load(data2, flag);
	Initial();
	Sort();
	vector<vector<float> >score = DisScore(m);
	vector<vector<int> >adjacent(NumBGene, vector<int>(NumBGene, 0));
	multimap<float, pair<int, int> >ma;
	for (int i = 0; i < NumBGene - 1; i++) {
		for (int j = i + 1; j < NumBGene; j++) {
			pair<int, int>te = make_pair(i, j);
			ma.insert(make_pair(score[i][j], te));
		}
	}
	int co = 0;
	int CIn = (NumBGene * (NumBGene - 1) / 2) * 0.001;
	multimap<float, pair<int, int> >::reverse_iterator iter;
	for (iter = ma.rbegin(); iter != ma.rend(); iter++) {
		if (co >= CIn) {
			break;
		}
		adjacent[iter->second.first][iter->second.second] = 1;
		adjacent[iter->second.second][iter->second.first] = 1;
		co++;
	}
	//DFS
	vector<int>visited(NumBGene, -1);
	int f = 0;
	for (int i = 0; i < NumBGene; i++) {
		if (visited[i] == -1) {
			visit(visited, i, f);
			DFS(adjacent, visited, i, f);
			f++;
		}
	}
	vector<vector<int> >count(f, vector<int>());
	for (int i = 0; i < NumBGene; i++) {
		count[visited[i]].push_back(i);
	}

	//generate new cluster
	for (int i = 0; i < count.size(); i++) {
		if (count[i].size() >= 3) {
			for (int j = 0; j < NumBGene; j++) {
				Group[j].push_back(0);
			}
			MarkedClusters.push_back(true);
			CountGene.push_back(count[i].size());
			Template.push_back(vector<int>());
			Record.push_back(vector<int>());
			TransferC.push_back(vector<float>());
			vector<float > temp(NumBCell, 0);
			for (int j = 0; j < NumBCell; j++) {
				Template[GeneClusters].push_back(0);
				Record[GeneClusters].push_back(0);
				TransferC[GeneClusters].push_back(0);
			}
			for (int j = 0; j < count[i].size(); j++) {
				Group[count[i][j]][GeneClusters] = 1;
				ClusterIndex[count[i][j]] = GeneClusters;
				for (int p = 0; p < NumBCell; p++) {
					Record[GeneClusters][p] = Record[GeneClusters][p] + BData[count[i][j]][p];
					temp[p] = temp[p] + Transfer[count[i][j]][p];
				}
			}
			for (int j = 0; j < NumBCell; j++) {
				if ((CountGene[GeneClusters] - Record[GeneClusters][j]) > Record[GeneClusters][j]) {//2 times
					Template[GeneClusters][j] = 0;
				}
				else {
					Template[GeneClusters][j] = 1;
				}
				TransferC[GeneClusters][j] = temp[j] / Record[GeneClusters][j];
			}
			GeneClusters++;
		}
	}
	for (int i = 0; i < MaxIterate; i++) {
		WeightedRho(i);
		for (int j = 0; j < NumBGene; j++) {
			RunGene(j);
		}
		Prune();
	}

	vector<vector<int> >temi;
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			temi.push_back(Template[i]);
		}
	}
	Second();
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			BLM.push_back(Template[i]);
		}
	}
	vector<vector<int> >temp(BLM[0].size(), vector<int>());
	for (int i = 0; i < BLM[0].size(); i++) {
		for (int j = 0; j < BLM.size(); j++) {
			temp[i].push_back(BLM[j][i]);
		}
	}
	
	return temp;
}

void	BG::Second() {
	SC();
	vector<vector<int> >temi;
	vector<vector<float> >te;
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			temi.push_back(Template[i]);
			te.push_back(TransferC[i]);
		}
	}
	for (int i = 0; i < GeneClusters; i++) {
		Template.clear();
		CountGene.clear();
		Record.clear();
		TransferC.clear();
	}
	MarkedClusters.clear();
	Template = temi;
	TransferC = te;
	for (int i = 0; i < Template.size(); i++) {
		MarkedClusters.push_back(true);
		CountGene.push_back(0);
		Record.push_back(vector<int>());
		for (int j = 0; j < NumBCell; j++) {
			Record[i].push_back(Template[i][j]);
		}
	}
	GeneClusters = Template.size();
	for (int i = 0; i < NumBGene; i++) {
		Group[i].clear();
		for (int j = 0; j < GeneClusters; j++) {
			Group[i].push_back(0);
		}
	}
	for (int i = 0; i < NumBGene; i++) {
		ClusterIndex[i] = -1;
	}
	for (int i = 0; i < MaxIterate; i++) {
		WeightedRho(i);
		for (int j = 0; j < NumBGene; j++) {
			X = BData[j];
			int p = RunComS(j);
			TrainS(j, p);
			X.clear();
		}
		Prune();
	}
}

void BG::visit(vector<int>& visited, int num, int flag) {
	visited[num] = flag;
}

void BG::DFS(vector<vector<int> >& adj, vector<int>& visited, int num, int flag) {
	for (int i = num + 1; i < NumBGene; i++) {
		if (adj[num][i] == 1) {
			if (visited[i] == -1) {
				visit(visited, i, flag);
				DFS(adj, visited, i, flag);
			}
		}
	}
}

void BG::SC() {
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			Sign.push_back(vector<float>(NumBCell, 0));
			for (int j = 0; j < NumBCell; j++) {
				if (Template[i][j] == 1) {
					Sign[Sign.size() - 1][j] = fabs(Record[i][j] / CountGene[i] - Per[j]) / Per[j];
					Sign[Sign.size() - 1][j] = 1 / (1 + exp(Sign[Sign.size() - 1][j]));
				}
				else {
					Sign[Sign.size() - 1][j] = fabs((CountGene[i] - Record[i][j]) / CountGene[i] - 1 + Per[j]) / (1 - Per[j]);
					Sign[Sign.size() - 1][j] = 1 / (1 + exp(Sign[Sign.size() - 1][j]));
				}
			}
		}
	}
}


