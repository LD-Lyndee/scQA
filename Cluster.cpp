
#include"Cluster.h"

CC::CC() {

}

void CC::Load(vector<vector<int> >& data1, vector<vector<float> >& data2) {
	for (int i = 0; i < data1.size(); i++) {
		BLM.push_back(data1[i]);
	}
	for (int i = 0; i < data2.size(); i++) {
		QLM.push_back(data2[i]);
	}
	NumCell = BLM.size();
	NumB = BLM[0].size();
	NumQ = QLM[0].size();
}

vector<vector<float> > CC::SimScore(char& m, int n) {
	vector<vector<float> >score(NumCell, vector<float>(NumCell, 0));
	vector<vector<float> >rank;
	switch (m) {
	case'M':
		//score with Manhattan Distance
		for (int i = 0; i < NumCell - 1; i++) {
			for (int j = i + 1; j < NumCell; j++) {
				if (n == 0) {
					for (int p = 0; p < NumB; p++) {
						score[i][j] = score[i][j] + fabs(BLM[i][p] - BLM[j][p]);
					}
					score[i][j] = 1 - score[i][j] / BLM[0].size();
				}
				else {
					for (int p = 0; p < NumQ; p++) {
						score[i][j] = score[i][j] + fabs(QLM[i][p] - QLM[j][p]);
					}
					score[i][j] = 1 - score[i][j] / QLM[0].size();
				}
				score[j][i] = score[i][j];
			}
		}
		break;
	case'C':
		//score with chi-square
		for (int i = 0; i < NumCell - 1; i++) {
			for (int j = i + 1; j < NumCell; j++) {
				int a = 0, b = 0, c = 0, d = 0;
				if (n == 0) {
					for (int p = 0; p < NumB; p++) {
						if (BLM[i][p] == 1 && BLM[j][p] == 1) {
							a++;
						}
						else if (BLM[i][p] == 0 && BLM[j][p] == 1) {
							b++;
						}
						else if (BLM[i][p] == 1 && BLM[j][p] == 0) {
							c++;
						}
						else
							d++;
					}
				}
				else {
					for (int p = 0; p < NumQ; p++) {
						if (QLM[i][p] >= 0.2 && QLM[j][p] >= 0.2) {
							a++;
						}
						else if (QLM[i][p] < 0.2 && QLM[j][p] >= 0.2) {
							b++;
						}
						else if (QLM[i][p] >= 0.2 && QLM[j][p] < 0.2) {
							c++;
						}
						else
							d++;
					}
				}
				score[i][j] = (float)((a * d - b * c) * (a * d - b * c)) / (float)((a + b) * (a + c) * (b + d) * (c + d));
				score[j][i] = score[i][j];
			}
		}
		break;
	case'R':
		//score with rand-index
		for (int i = 0; i < NumCell - 1; i++) {
			for (int j = i + 1; j < NumCell; j++) {
				int a = 0, b = 0;
				if (n == 0) {
					for (int p = 0; p < NumB; p++) {
						if ((BLM[i][p] == 1 && BLM[j][p] == 1) || (BLM[i][p] == 0 && BLM[j][p] == 0)) {
							a++;
						}
						else
							b++;
					}
				}
				else {
					for (int p = 0; p < NumQ; p++) {
						if ((QLM[i][p] >= 0.2 && QLM[j][p] >= 0.2) || (QLM[i][p] < 0.2 && QLM[j][p] < 0.2)) {
							a++;
						}
						else
							b++;
					}
				}
				score[i][j] = (float)a / (a + b);
				score[j][i] = score[i][j];
			}
		}
		break;
	case'F':
		//score with F1-score
		for (int i = 0; i < NumCell - 1; i++) {
			for (int j = i + 1; j < NumCell; j++) {
				int a = 0, b = 0, c = 0, d = 0;
				if (n == 0) {
					for (int p = 0; p < NumB; p++) {
						if (BLM[i][p] == 1 && BLM[j][p] == 1) {
							a++;
						}
						else if (BLM[i][p] == 0 && BLM[j][p] == 1) {
							b++;
						}
						else if (BLM[i][p] == 1 && BLM[j][p] == 0) {
							c++;
						}
						else
							d++;
					}
				}
				else {
					for (int p = 0; p < NumQ; p++) {
						if (QLM[i][p] >= 0.2 && QLM[j][p] >= 0.2) {
							a++;
						}
						else if (QLM[i][p] < 0.2 && QLM[j][p] >= 0.2) {
							b++;
						}
						else if (QLM[i][p] >= 0.2 && QLM[j][p] < 0.2) {
							c++;
						}
						else
							d++;
					}
				}
				score[i][j] = (float)(2 * a) / (float)(2 * a + b + c);
				score[j][i] = score[i][j];
			}
		}
		break;
	case'S':
		//score with spearman rank-score
		//sort, calculate rank
		if (n == 0) {
			for (int i = 0; i < NumCell; i++) {
				rank.push_back(vector<float>());
				vector<int>zero;
				vector<float>num(NumB, 0);
				for (int j = 0; j < NumB; j++) {
					if (BLM[i][j] == 0) {
						zero.push_back(j);
					}
				}
				//a1*n+n*(n-1)*d/2
				float sum = 1 + (zero.size() - 1) / 2;
				for (int j = 0; j < zero.size(); j++) {
					num[zero[j]] = sum;
				}
				sum = (zero.size() + 1) + (NumB - zero.size() - 1) / 2;
				for (int j = 0; j < NumB; j++) {
					if (num[j] == 0) {
						num[j] = sum;
					}
				}
				for (int j = 0; j < NumB; j++) {
					rank[i].push_back(num[j]);
				}
			}
			//calculate similarity
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					float d = 0;
					for (int p = 0; p < NumB; p++) {
						d = d + (rank[i][p] - rank[j][p]) * (rank[i][p] - rank[j][p]);
					}
					score[i][j] = 1 - 6 * d / (NumB * (NumB * NumB - 1));
					score[j][i] = score[i][j];
				}
			}
		}
		else {
			for (int i = 0; i < NumCell; i++) {
				rank.push_back(vector<float>());
				vector<int>zero;
				vector<float>num(NumQ, 0);
				multimap<float, int> map;
				for (int j = 0; j < NumQ; j++) {
					if (QLM[i][j] == 0) {
						zero.push_back(j);
					}
					else {
						map.insert(make_pair(QLM[i][j], j));
					}
				}
				//a1*n+n*(n-1)*d/2
				float sum = 1 + (zero.size() - 1) / 2;
				for (int j = 0; j < zero.size(); j++) {
					num[zero[j]] = sum;
				}
				sum = zero.size() + 1;
				multimap<float, int>::iterator iter;
				for (iter = map.begin(); iter != map.end(); iter++) {
					num[iter->second] = sum;
					sum++;
				}
				for (int j = 0; j < NumQ; j++) {
					rank[i].push_back(num[j]);
				}
			}
			//calculate similarity
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					float d = 0;
					for (int p = 0; p < NumQ; p++) {
						d = d + (rank[i][p] - rank[j][p]) * (rank[i][p] - rank[j][p]);
					}
					score[i][j] = 1 - 6 * d / (NumQ * (NumQ * NumQ - 1));
					score[j][i] = score[i][j];
				}
			}
		}
		break;
	case'P':
		//score with correlation 
		if (n == 0) {
			vector<float>scores(NumCell, 0);
			for (int i = 0; i < NumCell; i++) {
				for (int j = 0; j < NumB; j++) {
					scores[i] = scores[i] + BLM[i][j];
				}
				scores[i] = scores[i] / NumB;
			}
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					float sm = 0;
					float sn1 = 0;
					float sn2 = 0;
					for (int p = 0; p < NumB; p++) {
						sm = sm + (BLM[i][p] - scores[i]) * (BLM[j][p] - scores[j]);
						sn1 = sn1 + (BLM[i][p] - scores[i]) * (BLM[i][p] - scores[i]);
						sn2 = sn2 + (BLM[j][p] - scores[j]) * (BLM[j][p] - scores[j]);
					}
					score[i][j] = sm / (sqrt(sn1 * sn2));
					score[j][i] = score[i][j];
				}
			}
		}
		else {
			vector<float>scores(NumCell, 0);
			for (int i = 0; i < NumCell; i++) {
				for (int j = 0; j < NumQ; j++) {
					scores[i] = scores[i] + QLM[i][j];
				}
				scores[i] = scores[i] / NumQ;
			}
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					float sm = 0;
					float sn1 = 0;
					float sn2 = 0;
					for (int p = 0; p < NumQ; p++) {
						sm = sm + (QLM[i][p] - scores[i]) * (QLM[j][p] - scores[j]);
						sn1 = sn1 + (QLM[i][p] - scores[i]) * (QLM[i][p] - scores[i]);
						sn2 = sn2 + (QLM[j][p] - scores[j]) * (QLM[j][p] - scores[j]);
					}
					score[i][j] = sm / (sqrt(sn1 * sn2));
					score[j][i] = score[i][j];
				}
			}
		}
		break;
	case'K':
		//score with KL
		if (n == 0) {
			vector<vector<float> >BLM2(NumCell, vector<float>(NumCell, 0));
			for (int i = 0; i < NumCell; i++) {
				float sum = 0;
				int c = 0;
				float ep = 0.001;
				for (int j = 0; j < NumB; j++) {
					if (BLM[i][j] == 0) {
						c++;
					}
					sum = sum + BLM[i][j];
				}
				for (int j = 0; j < NumB; j++) {
					if (BLM[i][j] == 0) {
						BLM2[i][j] = ep / c;
					}
					else {
						BLM2[i][j] = BLM[i][j] / sum - ep / (NumB - c);
					}
				}
			}
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					for (int p = 0; p < NumB; p++) {
						score[i][j] = score[i][j] + BLM2[i][p] * log(BLM2[i][p] / BLM2[j][p]);
						score[j][i] = score[j][i] + BLM2[j][p] * log(BLM2[j][p] / BLM2[i][p]);
					}
				}
			}
			for (int i = 0; i < NumCell; i++) {
				float sum = 0;
				for (int j = 0; j < NumCell; j++) {
					sum = sum + score[i][j];
				}
				for (int j = 0; j < NumCell; j++) {
					if (i != j) {
						score[i][j] = (sum - score[i][j]) / sum;
					}
				}
			}
			//¶Ô³Æ
			for (int i = 0; i < NumCell; i++) {
				for (int j = 0; j < NumCell; j++) {
					if (i != j) {
						score[i][j] = (score[i][j] + score[j][i]) / 2;
						score[j][i] = score[i][j];
					}
				}
			}
		}
		else {
			vector<vector<float> >QLM2(NumCell, vector<float>(NumCell, 0));
			for (int i = 0; i < NumCell; i++) {
				float sum = 0;
				int c = 0;
				float ep = 0.001;
				for (int j = 0; j < NumQ; j++) {
					if (QLM[i][j] == 0) {
						c++;
					}
					sum = sum + QLM[i][j];
				}
				for (int j = 0; j < NumQ; j++) {
					if (QLM[i][j] == 0) {
						QLM2[i][j] = ep / c;
					}
					else {
						QLM2[i][j] = QLM[i][j] / sum - ep / (NumQ - c);
					}
				}
			}
			for (int i = 0; i < NumCell - 1; i++) {
				for (int j = i + 1; j < NumCell; j++) {
					for (int p = 0; p < NumQ; p++) {
						score[i][j] = score[i][j] + QLM2[i][p] * log(QLM2[i][p] / QLM2[j][p]);
						score[j][i] = score[j][i] + QLM2[j][p] * log(QLM2[j][p] / QLM2[i][p]);
					}
				}
			}
			for (int i = 0; i < NumCell; i++) {
				float sum = 0;
				for (int j = 0; j < NumCell; j++) {
					sum = sum + score[i][j];
				}
				for (int j = 0; j < NumCell; j++) {
					if (i != j) {
						score[i][j] = (sum - score[i][j]) / sum;
					}
				}
			}
			//¶Ô³Æ
			for (int i = 0; i < NumCell; i++) {
				for (int j = 0; j < NumCell; j++) {
					if (i != j) {
						score[i][j] = (score[i][j] + score[j][i]) / 2;
						score[j][i] = score[i][j];
					}
				}
			}
		}
		break;
	}
	return score;
}

void CC::SimGraph() {
	//calculate similarity
	char m1 = 'M';
	char m2 = 'P';
	vector<vector<float> >score1 = SimScore(m1, 0);
	
	vector<vector<float> >score2 = SimScore(m2, 1);
	
	//construct graph	
	//gap define neighbor
	for (int i = 0; i < NumCell; i++) {
		G2.push_back(vector<float>());//adjacency matrix using N2
		g2.push_back(vector<float>());//adjacency matrix using n2
		N1.push_back(vector<int>());//BLM neighbors 
		N2.push_back(vector<int>());//QLM neighbors 
		for (int j = 0; j < NumCell; j++) {
			G2[i].push_back(0);
			g2[i].push_back(0);
		}
	}
	
	//select neighbors 
	for (int i = 0; i < NumCell; i++) {
		vector<float>temp = score1[i];
		sort(temp.begin(), temp.end());
		reverse(temp.begin(), temp.end());
		int  sec = 0;
		for (int j = 0; j < NumCell; j++) {
			if (temp[j] == temp[0]) {
				continue;
			}
			else {
				sec = j;
				break;
			}
		}
		for (int j = 0; j < score1[i].size(); j++) {
			if (score1[i][j] == temp[0]) {
				N1[i].push_back(j);
			}
		}
		temp = score2[i];
		sort(temp.begin(), temp.end());
		reverse(temp.begin(), temp.end());
		
		int count;
		if (NumCell < 10000) {
			count = 0.1 * NumCell;
		}
		else {
			count = 0.05 * NumCell;
		}
		
		for (int j = 0; j < score2[i].size(); j++) {
			if (score2[i][j] >= temp[count]) {
				N2[i].push_back(j);
				G2[i][j] = score2[i][j];
			}
			else if (score2[i][j] >= temp[2 * count]) {
				g2[i][j] = score2[i][j];
			}
		}

	}
	PA();
}

void CC::PA() {
	int t = 2;
	vector<vector<int> >BGroup;
	vector<int>Label(NumCell, -1);
	vector<int>fb(NumCell, -1);
	if (NumCell < 500) {
		NumB = BLM[0].size();
		int MaxIterator = 6;
		vector<vector<int> >BCluster;
		vector<int>BF;
		vector<vector<int> >Record;
		for (int i = 0; i < MaxIterator; i++) {
			float Rho = 0.9 - 0.05 * i;
			for (int j = 0; j < NumCell; j++) {
				if (fb[j] == -1) {
					if (BCluster.size() == 0) {
						BCluster.push_back(BLM[j]);
						BF.push_back(1);
						Record.push_back(BLM[j]);
						fb[j] = 0;
					}
					else {
						multimap<float, int>CluSim;
						for (int p = 0; p < BCluster.size(); p++) {
							if (BF[p] >= 1) {
								float sim = 0;
								for (int q = 0; q < NumB; q++) {
									sim = sim + fabs(BLM[j][q] - BCluster[p][q]);
								}
								sim = sim / NumB;
								CluSim.insert(make_pair(sim, p));
							}
						}
						multimap<float, int>::iterator iter = CluSim.begin();
						if (CluSim.size() == 0) {
							BCluster.push_back(BLM[j]);
							Record.push_back(BLM[j]);
							BF.push_back(1);
							fb[j] = BF.size() - 1;
						}
						else {
							if (1 - iter->first < Rho) {
								BCluster.push_back(BLM[j]);
								Record.push_back(BLM[j]);
								BF.push_back(1);
								fb[j] = BF.size() - 1;
							}
							else {
								for (int p = 0; p < NumB; p++) {
									if (BLM[j][p] == 1) {
										Record[iter->second][p]++;
									}
									if (Record[iter->second][p] >= BF[iter->second] - Record[iter->second][p]) {
										BCluster[iter->second][p] = 1;
									}
									else {
										BCluster[iter->second][p] = 0;
									}
								}
								BF[iter->second]++;
								fb[j] = iter->second;
							}
						}
					}
				}
			}
			if (i != MaxIterator - 1) {
				for (int j = 0; j < BF.size(); j++) {
					if (BF[j] < NumCell * 0.1 && BF[j] != -1) {
						BF[j] = -1;
						for (int p = 0; p < NumCell; p++) {
							if (fb[p] == j) {
								fb[p] = -1;
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < BF.size(); i++) {
			if (BF[i] != -1) {
				BGroup.push_back(vector<int>());
				for (int j = 0; j < NumCell; j++) {
					if (fb[j] == i) {
						BGroup[BGroup.size() - 1].push_back(j);
					}
				}
			}
		}
		for (int i = 0; i < NumCell; i++) {
			fb[i] = -1;
		}
		for (int i = 0; i < BGroup.size(); i++) {
			for (int j = 0; j < BGroup[i].size(); j++) {
				fb[BGroup[i][j]] = i;
			}
		}
		N1.clear();
		for (int i = 0; i < NumCell; i++) {
			N1.push_back(vector<int>());
			if (fb[i] == -1)
			{
				continue;
			}
			for (int j = 0; j < BGroup[fb[i]].size(); j++) {
				if (i != BGroup[fb[i]][j]) {
					N1[i].push_back(BGroup[fb[i]][j]);
				}
			}
		}
	}
	else {
		for (int i = 0; i < NumCell; i++) {
			if (fb[i] != -1) {
				continue;
			}
			BGroup.push_back(vector<int>());
			BGroup[BGroup.size() - 1].push_back(i);
			fb[i] = BGroup.size() - 1;
			for (int j = 0; j < N1[i].size(); j++) {
				if (fb[N1[i][j]] == -1) {
					BGroup[BGroup.size() - 1].push_back(N1[i][j]);
					fb[N1[i][j]] = BGroup.size() - 1;
				}
			}
		}
		if (NumCell > 20000) {
			t = 1;
		}
	}

	vector<vector<float> >SimQ(NumCell, vector<float>(NumCell, 0));
	vector<int>BDeg;
	vector<vector<int> >BN(NumCell, vector<int>());
	vector<vector<int> >RN(NumCell, vector<int>());

	for (int i = 0; i < NumCell; i++) {
		for (int j = 0; j < N2[i].size(); j++) {
			RN[N2[i][j]].push_back(i);
		}
	}

	for (int i = 0; i < NumCell; i++) {
		for (int j = 0; j < NumCell; j++) {
			if (G2[i][j] != 0 && G2[j][i] != 0) {
				BN[i].push_back(j);
			}
		}
		BDeg.push_back(BN[i].size());
	}

	vector<float> SD(NumCell, 0);
	vector<multimap<float, int> >SD2(BGroup.size());
	multimap<float, int > SD3;
	//seed determination
	//compute average similarity
	for (int i = 0; i < NumCell; i++) {
		if (N2[i].size() == 0) {
			SD[i] = 0;
			SD2[fb[i]].insert(make_pair(0, i));
			SD3.insert(make_pair(0, i));
			continue;
		}
		float sum = 0;
		for (int j = 0; j < N2[i].size(); j++) {
			sum = sum + G2[i][N2[i][j]];
		}
		SD[i] = sum / N2[i].size();
		SD2[fb[i]].insert(make_pair(SD[i], i));
		SD3.insert(make_pair(SD[i], i));
	}

	//every group must have seed
	vector<float>mseed;
	for (int i = 0; i < BGroup.size(); i++) {
		multimap<float, int>::reverse_iterator iter = SD2[i].rbegin();
		mseed.push_back(iter->first);
	}
	sort(mseed.begin(), mseed.end());

	vector<vector<int> >CCS;
	if (NumCell < 500) {
		for (int i = 0; i < BGroup.size(); i++) {
			if (BGroup[i].size() <= 2) {
				continue;
			}
			vector<int>GSD;
			vector<int>Lab;
			for (int j = 0; j < BGroup[i].size(); j++) {
				if (SD[BGroup[i][j]] >= mseed[0]) {
					GSD.push_back(BGroup[i][j]);
					Lab.push_back(GSD.size());
				}
			}
			//one order neighbors
			if (GSD.size() > 1) {
				for (int j = 1; j < GSD.size(); j++) {
					for (int p = 0; p < j; p++) {
						if (G2[GSD[j]][GSD[p]] > 0 && G2[GSD[p]][GSD[j]] > 0) {
							Lab[j] = Lab[p];
						}
					}
				}
				vector<vector<int> >tCCS;
				tCCS.push_back(vector<int>());
				tCCS[tCCS.size() - 1].push_back(0);
				for (int j = 1; j < GSD.size(); j++) {
					int a = 0;
					for (int p = 0; p < tCCS.size(); p++) {
						if (Lab[j] == Lab[tCCS[p][0]]) {
							tCCS[p].push_back(j);
							a = 1;
							break;
						}
					}
					if (a == 0) {
						tCCS.push_back(vector<int>());
						tCCS[tCCS.size() - 1].push_back(j);
					}
				}
				for (int j = 0; j < tCCS.size(); j++) {
					CCS.push_back(vector<int>());
					for (int p = 0; p < tCCS[j].size(); p++) {
						CCS[CCS.size() - 1].push_back(GSD[tCCS[j][p]]);
					}
				}
			}
			else {
				CCS.push_back(vector<int>());
				CCS[CCS.size() - 1].push_back(GSD[0]);
			}
		}
	}
	else {
		for (int i = 0; i < BGroup.size(); i++) {
			if (BGroup[i].size() < 0.01 * NumCell) {
				continue;
			}
			vector<int>GSD;
			vector<int>Lab;
			for (int j = 0; j < BGroup[i].size(); j++) {
				if (SD[BGroup[i][j]] >= mseed[0]) {
					GSD.push_back(BGroup[i][j]);
					Lab.push_back(GSD.size());
				}
			}
			//one order neighbors
			if (GSD.size() > 1) {
				for (int j = 1; j < GSD.size(); j++) {
					for (int p = 0; p < j; p++) {
						if (G2[GSD[j]][GSD[p]] > 0 && G2[GSD[p]][GSD[j]] > 0) {
							Lab[j] = Lab[p];
						}
					}
				}
				vector<vector<int> >tCCS;
				tCCS.push_back(vector<int>());
				tCCS[tCCS.size() - 1].push_back(0);
				for (int j = 1; j < GSD.size(); j++) {
					int a = 0;
					for (int p = 0; p < tCCS.size(); p++) {
						if (Lab[j] == Lab[tCCS[p][0]]) {
							tCCS[p].push_back(j);
							a = 1;
							break;
						}
					}
					if (a == 0) {
						tCCS.push_back(vector<int>());
						tCCS[tCCS.size() - 1].push_back(j);
					}
				}
				for (int j = 0; j < tCCS.size(); j++) {
					CCS.push_back(vector<int>());
					for (int p = 0; p < tCCS[j].size(); p++) {
						CCS[CCS.size() - 1].push_back(GSD[tCCS[j][p]]);
					}
				}
			}
		}
	}
	//assign other cells
	float alpha = 0.5;
	for (int i = 0; i < CCS.size(); i++) {
		for (int j = 0; j < CCS[i].size(); j++) {
			Label[CCS[i][j]] = i;
		}
	}

	//change
	while (t > 0) {
		t--;
		for (multimap<float, int>::reverse_iterator iter = SD3.rbegin(); iter != SD3.rend(); iter++) {
			if (Label[iter->second] != -1) {
				continue;
			}
			int tp = iter->second;
			vector<vector<float> >temp(3, vector<float>(CCS.size(), 0));
			vector<vector<int> >nt(2, vector<int>(CCS.size(), 0));
			float unlabel = 0;
			int nun = 0;
			float runlabel = 0;
			int nrun = 0;
			for (int j = 0; j < N2[tp].size(); j++) {
				if (Label[N2[tp][j]] != -1) {
					temp[0][Label[N2[tp][j]]] = temp[0][Label[N2[tp][j]]] + G2[tp][N2[tp][j]];
					nt[0][Label[N2[tp][j]]]++;
				}
				else {
					unlabel = unlabel + G2[tp][N2[tp][j]];
					nun++;
					for (int p = 0; p < N2[N2[tp][j]].size(); p++) {
						if (Label[N2[N2[tp][j]][p]] != -1) {
							temp[0][Label[N2[N2[tp][j]][p]]] = temp[0][Label[N2[N2[tp][j]][p]]] + G2[N2[tp][j]][N2[N2[tp][j]][p]];
							nt[0][Label[N2[N2[tp][j]][p]]]++;
						}
					}
				}
			}
			for (int j = 0; j < RN[tp].size(); j++) {
				if (Label[RN[tp][j]] != -1) {
					temp[1][Label[RN[tp][j]]] = temp[1][Label[RN[tp][j]]] + G2[tp][RN[tp][j]];
					nt[1][Label[RN[tp][j]]]++;
				}
				else {
					runlabel = runlabel + G2[tp][RN[tp][j]];
					nrun++;
					for (int p = 0; p < RN[RN[tp][j]].size(); p++) {
						if (Label[RN[RN[tp][j]][p]] != -1) {
							temp[0][Label[RN[RN[tp][j]][p]]] = temp[0][Label[RN[RN[tp][j]][p]]] + G2[RN[tp][j]][RN[RN[tp][j]][p]];
							nt[0][Label[RN[RN[tp][j]][p]]]++;
						}
					}
				}
			}
			multimap<float, int>st;
			for (int j = 0; j < CCS.size(); j++) {
				if (nt[0][j] == 0) {
					temp[2][j] = 0;
				}
				else {
					temp[2][j] = alpha * temp[0][j] / nt[0][j];
				}
				if (nt[1][j] == 0) {
					temp[2][j] = temp[2][j] + 0;
				}
				else {
					temp[2][j] = temp[2][j] + (1 - alpha) * temp[1][j] / nt[1][j];
				}
				st.insert(make_pair(temp[2][j], j));
			}
			if (NumCell < 500) {
				float u = 0;
				if (nun == 0) {
					if (nrun == 0) {
						u = 0;
					}
					else {
						u = u + (1 - alpha) * runlabel / nrun;
					}
				}
				else {
					if (nrun == 0) {
						u = alpha * unlabel / nun;
					}
					else {
						u = alpha * unlabel / nun + (1 - alpha) * runlabel / nrun;
					}
				}
				multimap<float, int>::reverse_iterator it = st.rbegin();
				if (it->first > u) {
					CCS[it->second].push_back(tp);
				}
			}
			else {
				if (st.size() > 0) {
					multimap<float, int>::reverse_iterator it = st.rbegin();
					if (it->first > alpha * unlabel / nun + (1 - alpha) * runlabel / nrun) {
						CCS[it->second].push_back(tp);
					}
				}
			}

		}
	}
	//delete community of size less than numcell*0.01
	vector<vector<int> >tccs;
	if (NumCell > 20000) {
		vector<vector<int> >tccs;
		for (int i = 0; i < CCS.size(); i++) {
			if (CCS[i][0] != -1) {
				tccs.push_back(CCS[i]);
			}
		}
		CCS = tccs;
		for (int i = 0; i < NumCell; i++) {
			Label[i] = -1;
		}
		for (int i = 0; i < CCS.size(); i++) {
			for (int j = 0; j < CCS[i].size(); j++) {
				Label[CCS[i][j]] = i;
			}
		}
		vector<int>unlabel;
		for (int i = 0; i < NumCell; i++) {
			if (Label[i] == -1) {
				unlabel.push_back(i);
			}
		}
		vector<float>ave;
		for (int i = 0; i < CCS.size(); i++) {
			vector<int>nei(CCS[i].size(), 0);
			for (int j = 0; j < CCS[i].size(); j++) {
				for (int p = 0; p < N2[CCS[i][j]].size(); p++) {
					if (Label[N2[CCS[i][j]][p]] == i) {
						nei[j]++;
					}
				}
			}
			float sum = 0;
			for (int j = 0; j < nei.size(); j++) {
				sum = sum + nei[j];
			}
			sum = sum / CCS[i].size();
			ave.push_back(sum);
		}

		multimap<int, int>untemp;
		for (int i = 0; i < unlabel.size(); i++) {
			vector<float>compare(CCS.size(), 0);
			for (int j = 0; j < N2[unlabel[i]].size(); j++) {
				if (Label[N2[unlabel[i]][j]] != -1) {
					compare[Label[N2[unlabel[i]][j]]]++;
				}
			}
			multimap<float, int>scom;
			for (int j = 0; j < CCS.size(); j++) {
				scom.insert(make_pair(compare[j] - ave[j], j));
			}
			multimap<float, int>::reverse_iterator iter = scom.rbegin();
			Label[unlabel[i]] = iter->second;
		}
	}
	else {
		float p1 = 0.01 * NumCell;
		if (NumCell < 500) {
			p1 = 0.1 * NumCell;
		}
		for (int i = 0; i < CCS.size(); i++) {
			if (CCS[i].size() < p1) {
				for (int j = 0; j < CCS[i].size(); j++) {
					Label[CCS[i][j]] = -1;
				}
				CCS[i].erase(CCS[i].begin(), CCS[i].end());
			}
			else {
				tccs.push_back(CCS[i]);
			}
		}

		CCS.clear();
		multimap<int, int>ssize;
		for (int i = 0; i < tccs.size(); i++) {
			ssize.insert(make_pair(tccs[i].size(), i));
		}
		for (multimap<int, int>::iterator iter = ssize.begin(); iter != ssize.end(); iter++) {
			CCS.push_back(tccs[iter->second]);
		}
		for (int i = 0; i < NumCell; i++) {
			Label[i] = -1;
		}
		for (int i = 0; i < CCS.size(); i++) {
			for (int j = 0; j < CCS[i].size(); j++) {
				Label[CCS[i][j]] = i;
			}
		}
		//merge community
		vector<vector<int> >MC(CCS.size(), vector<int>(CCS.size(), 0));
		for (int i = 0; i < CCS.size(); i++) {
			for (int j = 0; j < CCS[i].size(); j++) {
				for (int p = 0; p < N2[CCS[i][j]].size(); p++) {
					if (Label[N2[CCS[i][j]][p]] != -1) {
						MC[i][Label[N2[CCS[i][j]][p]]]++;
					}
				}
			}
		}
		vector<int>Merge;
		for (int i = 0; i < MC.size(); i++) {
			int max = 0;
			int id = 0;
			for (int j = 0; j < MC.size(); j++) {
				if (MC[i][j] > max) {
					max = MC[i][j];
					id = j;
				}
			}
			if (max > MC[i][i]) {
				Merge.push_back(id);
			}
			else {
				Merge.push_back(-1);
			}
		}

		vector<int>flag(CCS.size(), -1);
		for (int i = 0; i < CCS.size(); i++) {
			if (Merge[i] != -1 && flag[Merge[i]] == -1) {
				for (int j = 0; j < CCS[i].size(); j++) {
					CCS[Merge[i]].push_back(CCS[i][j]);
				}
				CCS[i][0] = -1;
				flag[i] = Merge[i];
				for (int j = 0; j < CCS.size(); j++) {
					if (flag[j] == i) {
						flag[j] = Merge[i];
					}
				}
			}
			else if (NumCell < 500&&Merge[i] != -1 && flag[Merge[i]] != -1 && flag[Merge[i]] != i) {
				for (int j = 0; j < CCS[i].size(); j++) {
					CCS[flag[Merge[i]]].push_back(CCS[i][j]);
				}
				CCS[i][0] = -1;
				flag[i] = flag[Merge[i]];
			}
			else if (NumCell >= 500 && Merge[i] != -1 && flag[Merge[i]] != -1 ) {
				for (int j = 0; j < CCS[i].size(); j++) {
					CCS[flag[Merge[i]]].push_back(CCS[i][j]);
				}
				CCS[i][0] = -1;
				flag[i] = flag[Merge[i]];
			}
		}
		for (int i = 0; i < CCS.size(); i++) {
			if (CCS[i][0] == -1) {
				continue;
			}
			for (int j = 0; j < CCS[i].size(); j++) {
				Label[CCS[i][j]] = i;
			}
		}
		//unlabel 
		//average neighbors
		tccs.clear();
		for (int i = 0; i < CCS.size(); i++) {
			if (CCS[i][0] != -1) {
				tccs.push_back(CCS[i]);
			}
		}
		CCS = tccs;
		for (int i = 0; i < NumCell; i++) {
			Label[i] = -1;
		}
		for (int i = 0; i < CCS.size(); i++) {
			for (int j = 0; j < CCS[i].size(); j++) {
				Label[CCS[i][j]] = i;
			}
		}
		vector<int>unlabel;
		for (int i = 0; i < NumCell; i++) {
			if (Label[i] == -1) {
				unlabel.push_back(i);
			}
		}
		vector<float>ave;
		for (int i = 0; i < CCS.size(); i++) {
			vector<int>nei(CCS[i].size(), 0);
			for (int j = 0; j < CCS[i].size(); j++) {
				for (int p = 0; p < N2[CCS[i][j]].size(); p++) {
					if (Label[N2[CCS[i][j]][p]] == i) {
						nei[j]++;
					}
				}
			}
			float sum = 0;
			for (int j = 0; j < nei.size(); j++) {
				sum = sum + nei[j];
			}
			sum = sum / CCS[i].size();
			ave.push_back(sum);
		}

		multimap<int, int>untemp;
		for (int i = 0; i < unlabel.size(); i++) {
			vector<float>compare(CCS.size(), 0);
			for (int j = 0; j < N2[unlabel[i]].size(); j++) {
				if (Label[N2[unlabel[i]][j]] != -1) {
					compare[Label[N2[unlabel[i]][j]]]++;
				}
			}
			multimap<float, int>scom;
			for (int j = 0; j < CCS.size(); j++) {
				scom.insert(make_pair(compare[j] - ave[j], j));
			}
			multimap<float, int>::reverse_iterator iter = scom.rbegin();
			if (iter->first > 0&& scom.size() > 0) {
				Label[unlabel[i]] = iter->second;
			}
			else {
				int sum = N2[unlabel[i]].size();
				for (int j = 0; j < CCS.size(); j++) {
					sum = sum - compare[j];
				}
				untemp.insert(make_pair(sum, unlabel[i]));
			}
		}
		//vector<int>L3 = Label;
		//unlabel
		//new cluster
		vector<vector<int> >New;
		for (multimap<int, int>::reverse_iterator iter = untemp.rbegin(); iter != untemp.rend(); iter++) {
			if (New.size() == 0) {
				New.push_back(vector<int>());
				New[0].push_back(iter->second);
				Label[iter->second] = CCS.size();
				int i = iter->second;
				for (int j = 0; j < N2[i].size(); j++) {
					if (G2[i][N2[i][j]] > 0 && G2[N2[i][j]][i] > 0 && Label[N2[i][j]] == -1) {
						New[0].push_back(N2[i][j]);
						Label[N2[i][j]] = CCS.size();
					}
				}
			}
		}

		if (New.size() != 0) {
			CCS.push_back(New[0]);
		}
		int d = 1;
		while (d) {
			int c = d;
			d = 0;
			for (int i = 0; i < NumCell; i++) {
				if (Label[i] != -1) {
					continue;
				}
				vector<int>nei(CCS.size(), 0);
				for (int j = 0; j < N2[i].size(); j++) {
					if (Label[N2[i][j]] != -1) {
						nei[Label[N2[i][j]]]++;
					}
				}
				int max = 0;
				int id = -1;
				for (int j = 0; j < CCS.size(); j++) {
					if (nei[j] > max) {
						max = nei[j];
						id = j;
					}
				}
				Label[i] = id;
				if (id == -1) {
					d = i;
				}
			}
			if (c == d&&NumCell>=500) {
				break;
			}
		}
		for (int i = 0; i < CCS.size(); i++) {
			CCS[i].clear();
		}
		for (int i = 0; i < NumCell; i++) {
			if (Label[i] != -1) {
				CCS[Label[i]].push_back(i);
			}
		}


		//prune
		//search wrong node
		
		vector<vector<int> >RCS(CCS.size(), vector<int>());
		int c = 0;
		vector<int>L2 = Label;
		while (c == 0) {
			c = 1;
			for (int i = 0; i < RCS.size(); i++) {
				RCS[i].clear();
			}
			for (int i = 0; i < NumCell; i++) {
				vector<float>temp(CCS.size(), 0);
				for (int j = 0; j < N2[i].size(); j++) {
					temp[L2[N2[i][j]]] = temp[L2[N2[i][j]]] + G2[i][N2[i][j]];
					
				}
				float max = 0;
				int id = 0;
				for (int j = 0; j < CCS.size(); j++) {
					if (temp[j] > max) {
						max = temp[j];
						id = j;
					}
				}
				if (id != L2[i]) {
					L2[i] = id;
					c = 0;
				}
				RCS[L2[i]].push_back(i);
			}
		}
		
		if (NumCell < 500) {
			//merge
			MC.clear();
			CCS.clear();
			for (int i = 0; i < NumCell; i++) {
				if (CCS.size() == 0) {
					CCS.push_back(vector<int>());
					CCS[0].push_back(i);
				}
				else {
					int e = 0;
					for (int j = 0; j < CCS.size(); j++) {
						if (Label[i] == Label[CCS[j][0]]) {
							CCS[j].push_back(i);
							e = 1;
							break;
						}
					}
					if (e == 0) {
						CCS.push_back(vector<int>());
						CCS[CCS.size() - 1].push_back(i);
					}
				}
			}
			for (int i = 0; i < CCS.size(); i++) {
				for (int j = 0; j < CCS[i].size(); j++) {
					Label[CCS[i][j]] = i;
				}
			}
			MC.clear();
			for (int i = 0; i < CCS.size(); i++) {
				MC.push_back(vector<int>());
				for (int j = 0; j < CCS.size(); j++) {
					MC[i].push_back(0);
				}
			}
			for (int i = 0; i < CCS.size(); i++) {
				for (int j = 0; j < CCS[i].size(); j++) {
					for (int p = 0; p < N2[CCS[i][j]].size(); p++) {
						if (Label[N2[CCS[i][j]][p]] != -1) {
							MC[i][Label[N2[CCS[i][j]][p]]]++;
						}
					}
				}
			}
			Merge.clear();
			for (int i = 0; i < MC.size(); i++) {
				int max = 0;
				int id = 0;
				for (int j = 0; j < MC.size(); j++) {
					if (MC[i][j] > max) {
						max = MC[i][j];
						id = j;
					}
				}
				if (max > MC[i][i]) {
					Merge.push_back(id);
				}
				else {
					Merge.push_back(-1);
				}
			}
		}
	}
	FILE* fp = fopen("labels.txt", "wt");
	if (!fp) {
		printf("The file was not opened\n");
		exit(1);//Òì³£ÍË³ö
	}
	for (int i = 0; i < NumCell; i++) {
		fprintf(fp, "%d", Label[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void CC::Run(vector<vector<int> >& data1, vector<vector<float> >& data2) {
	Load(data1, data2);
	SimGraph();
}
