//quantitative cluster of genes

#include"Quantitative.h"

QG::QG() {
	RESET = 1;
	Rho = 0.7;
	GeneClusters = 0;
	MaxIterate = 1;
	MinGenes = 3;//least number of genes in a cluster
}

//preprocessing matrix, normalize to 0-1
vector<vector<float> >QG::PreMat(vector<vector<float> >& data) {
	vector<vector<float> >temp;
	for (int i = 0; i < data.size(); i++) {
		temp.push_back(vector<float>());
	}
	for (int i = 0; i < data[0].size(); i++) {
		float max = 0;
		for (int j = 0; j < data.size(); j++) {
			if (data[j][i] > max) {
				max = data[j][i];
			}
		}
		for (int j = 0; j < data.size(); j++) {
			temp[j].push_back(data[j][i] / max);
		}
	}
	return temp;
}

//load genes
void QG::Load(vector<vector<float> >& data, vector<int>& flag, vector<vector<int> >& RecordState) {
	for (int i = 0; i < flag.size(); i++) {
		QData.push_back(vector<float>());
		for (int j = 0; j < data[0].size(); j++) {
			QData[i].push_back(data[flag[i]][j]);
		}
	}
	NumQGene = QData.size();
	NumQCell = QData[0].size();
}

//initial vectors
void QG::Initial() {
	for (int i = 0; i < NumQGene; i++) {
		ClusterIndex.push_back(-1);
	}
	for (int i = 0; i < NumQGene; i++) {
		Group.push_back(vector<int>());
	}
}

void QG::WeightedRho(int j) {
	Rho = Rho - j * 0.05;
	Rho = max(float(0.7), Rho);
}

//initial clusters
vector<vector<float> > QG::SimilarityScore(char& m) {
	vector<vector<float> >score(NumQGene, vector<float>(NumQGene, 0));
	vector<float>scores(NumQGene, 0);
	vector<float>scores2(NumQGene, 0);
	vector<vector<float> >rank;
	vector<float> dd(NumQGene, 0);
	vector<vector<float> >QData2(NumQGene, vector<float>(NumQCell, 0));
	float ave = 0;
	float sm = 0;
	float sn1 = 0;
	float sn2 = 0;
	switch (m) {
	case'M':
		//score with Manhattan Distance
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				for (int p = 0; p < NumQCell; p++) {
					score[i][j] = score[i][j] + fabs(QData[i][p] - QData[j][p]);
				}
				score[i][j] = 1 - score[i][j] / NumQCell;
				score[j][i] = score[i][j];
			}
		}
		break;
	case'C':
		//score with chi-square
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				int a = 0, b = 0, c = 0, d = 0;
				for (int p = 0; p < NumQCell; p++) {
					if (QData[i][p] >= 0.2 && QData[j][p] >= 0.2) {
						a++;
					}
					else if (QData[i][p] < 0.2 && QData[j][p] >= 0.2) {
						b++;
					}
					else if (QData[i][p] >= 0.2 && QData[j][p] < 0.2) {
						c++;
					}
					else
						d++;
				}
				score[i][j] = (float)((a * d - b * c) * (a * d - b * c)) / (float)((a + b) * (a + c) * (b + d) * (c + d));
				score[j][i] = score[i][j];
			}
		}
		break;
	case'R':
		//score with rand-index
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				int a = 0, b = 0;
				for (int p = 0; p < NumQCell; p++) {
					if ((QData[i][p] >= 0.2 && QData[j][p] >= 0.2) || (QData[i][p] < 0.2 && QData[j][p] < 0.2)) {
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
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				int a = 0, b = 0, c = 0, d = 0;
				for (int p = 0; p < NumQCell; p++) {
					if (QData[i][p] >= 0.2 && QData[j][p] >= 0.2) {
						a++;
					}
					else if (QData[i][p] < 0.2 && QData[j][p] >= 0.2) {
						b++;
					}
					else if (QData[i][p] >= 0.2 && QData[j][p] < 0.2) {
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
	case'S':
		//score with spearman rank-score
		//sort, calculate rank
		for (int i = 0; i < NumQGene; i++) {
			rank.push_back(vector<float>());
			vector<int>zero;
			vector<float>num(NumQCell, 0);
			multimap<float, int> map;
			for (int j = 0; j < NumQCell; j++) {
				if (QData[i][j] == 0) {
					zero.push_back(j);
				}
				else {
					map.insert(make_pair(QData[i][j], j));
				}
			}
			float sum = 0;
			for (int j = 1; j < zero.size() + 1; j++) {
				sum = sum + j;
			}
			for (int j = 0; j < zero.size(); j++) {
				num[zero[j]] = sum / zero.size();
			}
			sum = zero.size() + 1;
			multimap<float, int>::iterator it;
			for (it = map.begin(); it != map.end(); it++) {
				num[it->second] = sum;
				sum++;
			}
			for (int j = 0; j < NumQCell; j++) {
				rank[i].push_back(num[j]);
			}
		}
		//calculate similarity
		ave = 0;
		for (int i = 1; i <= NumQCell; i++) {
			ave = ave + i;
		}
		ave = ave / NumQCell;
		for (int i = 0; i < NumQGene; i++) {
			for (int j = 0; j < NumQCell; j++) {
				dd[i] = dd[i] + (rank[i][j] - ave) * (rank[i][j] - ave);
			}
			dd[i] = sqrt(dd[i]);
		}
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				float m = 0;
				float d = dd[i] * dd[j];
				for (int p = 0; p < NumQCell; p++) {
					m = m + (rank[i][p] - ave) * (rank[j][p] - ave);
				}
				score[i][j] = m / d;
				score[j][i] = m / d;
			}
		}
		break;
	case'P':
		//score with correlation 
		for (int i = 0; i < NumQGene; i++) {
			for (int j = 0; j < NumQCell; j++) {
				scores[i] = scores[i] + QData[i][j];
			}
			scores[i] = scores[i] / NumQCell;
		}
		for (int i = 0; i < NumQGene - 1; i++) {
			for (int j = i + 1; j < NumQGene; j++) {
				sm = 0;//fenzi
				sn1 = 0;
				sn2 = 0;
				for (int p = 0; p < NumQCell; p++) {
					sm = sm + (QData[i][p] - scores[i]) * (QData[j][p] - scores[j]);
					sn1 = sn1 + (QData[i][p] - scores[i]) * (QData[i][p] - scores[i]);
					sn2 = sn2 + (QData[j][p] - scores[j]) * (QData[j][p] - scores[j]);
				}
				score[i][j] = sm / (sqrt(sn1 * sn2));
				score[j][i] = score[i][j];
			}
		}
		break;
	}
	return score;
}

void QG::Prune() {
	//Prune clusters
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			if (CountGene[i] < MinGenes) {
				MarkedClusters[i] = false;
				for (int j = 0; j < NumQGene; j++) {
					if (ClusterIndex[j] == i) {
						ClusterIndex[j] = -1;
						Group[j][i] = 0;
					}
				}
			}
		}
	}
}

void QG::visit(vector<int>& visited, int num, int flag) {
	visited[num] = flag;
}

void QG::DFS(vector<vector<int> >& adj, vector<int>& visited, int num, int flag) {
	for (int i = num + 1; i < NumQGene; i++) {
		if (adj[num][i] == 1) {
			if (visited[i] == -1) {
				visit(visited, i, flag);
				DFS(adj, visited, i, flag);
			}
		}
	}
}

vector<vector<float> > QG::QuantitativeCluster(vector<vector<float> >& data, vector<int>& flag, char& m, vector<vector<int> >& RecordState) {
	vector<vector<float> >temp = PreMat(data);
	Load(temp, flag, RecordState);
	Initial();
	vector<vector<float> >score = SimilarityScore(m);//calculate similarity
	multimap<float, pair<int, int> >ma;
	vector<vector<int> >adjacent(NumQGene, vector<int>(NumQGene, 0));
	for (int i = 0; i < NumQGene - 1; i++) {
		for (int j = i + 1; j < NumQGene; j++) {
			pair<int, int>te = make_pair(i, j);
			ma.insert(make_pair(score[i][j], te));
		}
	}
	int co = 0;
	multimap<float, pair<int, int> >::reverse_iterator iter;
	int CIn = (NumQGene * (NumQGene - 1) / 2) * 0.02;
	for (iter = ma.rbegin(); iter != ma.rend(); iter++) {
		if (co >= CIn) {
			break;
		}
		adjacent[iter->second.first][iter->second.second] = 1;
		adjacent[iter->second.second][iter->second.first] = 1;
		co++;
	}

	//DFS
	vector<int>visited(NumQGene, -1);
	int flag1 = 0;
	for (int i = 0; i < NumQGene; i++) {
		if (visited[i] == -1) {
			visit(visited, i, flag1);
			DFS(adjacent, visited, i, flag1);
			flag1++;
		}
	}
	vector<vector<int> >count(flag1, vector<int>());
	for (int i = 0; i < NumQGene; i++) {
		count[visited[i]].push_back(i);
	}
	//generate new cluster
	for (int i = 0; i < count.size(); i++) {
		if (count[i].size() > 2) {
			for (int j = 0; j < NumQGene; j++) {
				Group[j].push_back(0);
			}
			Ssort.push_back(vector<vector<float> >(NumQCell, vector<float>()));
			Template.push_back(vector<float>());
			Template2.push_back(vector<float>());
			Template3.push_back(vector<float>());
			for (int j = 0; j < NumQCell; j++) {
				Template[GeneClusters].push_back(0);
				Template2[GeneClusters].push_back(0);
				Template3[GeneClusters].push_back(1);
			}
			for (int j = 0; j < count[i].size(); j++) {
				Group[count[i][j]][GeneClusters] = 1;
				ClusterIndex[count[i][j]] = GeneClusters;
				for (int p = 0; p < NumQCell; p++) {
					Template[GeneClusters][p] = Template[GeneClusters][p] + QData[count[i][j]][p];
					Ssort[GeneClusters][p].push_back(QData[count[i][j]][p]);
				}
			}
			for (int j = 0; j < Ssort[GeneClusters].size(); j++) {
				sort(Ssort[GeneClusters][j].begin(), Ssort[GeneClusters][j].end());
				int low = 0.25 * Ssort[GeneClusters][j].size();
				int high = 0.75 * Ssort[GeneClusters][j].size();
				Template2[GeneClusters][j] = Ssort[GeneClusters][j][high];
				Template3[GeneClusters][j] = Ssort[GeneClusters][j][low];
			}
			for (int j = 0; j < NumQCell; j++) {
				Template[GeneClusters][j] = Template[GeneClusters][j] / count[i].size();
			}
			MarkedClusters.push_back(true);
			CountGene.push_back(count[i].size());
			GeneClusters++;
		}
	}

	for (int i = 0; i < MaxIterate; i++) {
		WeightedRho(i);
		for (int j = 0; j < NumQGene; j++) {
			if (ClusterIndex[j] == -1) {
				RunGene(j);
			}
		}
		Prune();
	}
	vector < vector<float> >temp2(NumQCell, vector<float>());
	for (int i = 0; i < GeneClusters; i++) {
		if (MarkedClusters[i]) {
			for (int j = 0; j < Template[i].size(); j++) {
				temp2[j].push_back(Template[i][j]);
			}
		}
	}
	return temp2;
}

void QG::RunGene(int tp) {
	X = QData[tp];
	int j = RunCom(tp);
	Train(tp, j);
	X.clear();
	RESET = 1;
}

int QG::RunCom(int tp) {
	if (GeneClusters == 0) {
		RESET = -1;
		return -1;
	}
	else {
		multimap<float, int>Choice;
		int index = 0;
		for (int i = 0; i < GeneClusters; i++) {
			if (MarkedClusters[i]) {
				float a = 0;
				for (int j = 0; j < NumQCell; j++) {
					if (X[j] > Template3[i][j] && X[j] < Template2[i][j]) {
						a = a + 1;
					}
				}
				Choice.insert(make_pair(a / NumQCell, i));
			}
		}
		while (RESET) {
			if (Choice.size() > 0) {
				multimap<float, int>::reverse_iterator iter = Choice.rbegin();
				if (iter->first >= Rho) {
					RESET = 0; //resonate
					return iter->second;
				}
				else {
					RESET = 1;
					return -1;
				}
			}
			else {
				RESET = -1;//new cluster
				return -1;
			}
		}
	}
}

void QG::Train(int tp, int j) {
	if (RESET == 0) {
		//resonate
		Group[tp][j] = 1;
		ClusterIndex[tp] = j;
		//update template
		for (int i = 0; i < NumQCell; i++) {
			Template[j][i] = (Template[j][i] * CountGene[j] + X[i]) / (CountGene[j] + 1);
			Ssort[j][i].push_back(X[i]);
			sort(Ssort[j][i].begin(), Ssort[j][i].end());
			int low = 0.25 * (CountGene[j] + 1);
			int high = 0.75 * (CountGene[j] + 1);
			Template2[j][i] = Ssort[j][i][high];
			Template3[j][i] = Ssort[j][i][low];
		}
		CountGene[j]++;
	}
	else {
		//RESET=-1, generate new cluster
		Generate(tp);
	}
}

void QG::Generate(int tp) {
	//generate new cluster
	Template.push_back(vector<float>());
	Template2.push_back(vector<float>());
	Template3.push_back(vector<float>());
	Ssort.push_back(vector<vector<float> >(NumQCell, vector<float>()));
	for (int i = 0; i < NumQCell; i++) {
		Template[GeneClusters].push_back(QData[tp][i]);
		Template2[GeneClusters].push_back(min(QData[tp][i] + 0.1, 1.0));
		Template3[GeneClusters].push_back(max(QData[tp][i] - 0.1, 0.0));
		Ssort[GeneClusters][i].push_back(QData[tp][i]);
	}
	for (int i = 0; i < NumQGene; i++) {
		Group[i].push_back(0);
	}
	Group[tp][GeneClusters] = 1;
	ClusterIndex[tp] = GeneClusters;
	MarkedClusters.push_back(true);
	CountGene.push_back(1);
	GeneClusters++;
}

