#include<iostream>
#include<math.h>
#include<string>
#include<limits.h>
#include<utility>
#include<vector>
#include<tuple>
#include<algorithm>
#include<stdlib.h>

using namespace std;
typedef long long ll;
			//          A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V        
int matrix[20][20] = { {4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0},
					  {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3},
				      {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1,0 ,-4,-2,-3},
				      {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3},
			       	  { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2 - 1},
				      {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2 ,1, 0,-3,-1, 0,-1,-2,-1,-2 },
				      {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2 },
					  { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3},
				      {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3},
				      {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3},
				      {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2 ,0,-3,-2,-1,-2,-1, 1},
				      {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2},
				      {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1,2 ,-1, 5, 0,-2,-1,-1,-1,-1 ,1},
				      {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0 ,0,-3, 0, 6,-4,-2,-2, 1 ,3,-1},
				      {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4 ,7,-1,-1,-4,-3,-2},
				      { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2} ,
				      { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0},
				      {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3},
				      {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1 ,3,-3,-2,-2, 2 ,7,-1 },
				      {0 ,-3,-3,-3,-1,-2,-2,-3,-3, 3 ,1,-2, 1,-1,-2,-2, 0,-3,-1, 4 } };

string protein = "ARNDCQEGHILKMFPSTWYV";
int arraynum(char a) {
	for (int i = 0; i < protein.size(); i++) {
		if (protein[i] == a) {
			return i;
		}
	}
	return -1;
}
tuple<int,string,string> score(string s1, string s2,int gap) {
	vector<vector<int>> SM(s1.size() + 1, vector<int>(s2.size() + 1, 0));
	vector<vector<string>>S1(s1.size() + 1, vector<string>(s2.size() + 1, s1));
	vector<vector<string>>S2(s1.size() + 1, vector<string>(s2.size() + 1, s2));
	vector<vector<int>>ins1(s1.size() + 1, vector<int>(s2.size() + 1));
	vector<vector<int>>ins2(s1.size() + 1, vector<int>(s2.size() + 1));
	ins1[0][0] = 0;
	ins2[0][0] = 0;
	S1[0][0] = s1;
	S2[0][0] = s2;
	for (int i = 1; i <= s1.size(); i++) {
		SM[i][0] = gap * i;
		S2[i][0] = S2[i - 1][0];
		S2[i][0].insert(0, "-");
		ins2[i][0] = ins2[i - 1][0] + 1;
		ins1[i][0] = ins1[i - 1][0] + 1;
	}
	for (int j = 1; j <= s2.size(); j++) {
		SM[0][j] = gap * j;
		S1[0][j] = S1[0][j - 1];
		S1[0][j].insert(0, "-");
		ins2[0][j] = ins2[0][j] + 1;
		ins1[0][j] = ins1[0][j] + 1;
	}	
	int a, b;
	for (int i = 1; i <= s1.size(); i++) {
		a = arraynum(s1[i - 1]);
		for (int j = 1; j <= s2.size(); j++) {
			b = arraynum(s2[j - 1]);
			int tmp = 0;
			SM[i][j] = SM[i][j - 1] + gap;
			if (a == -1) {
				SM[i][j] = SM[i][j - 1];
			}
			if (b == -1) {
				if (SM[i][j] < SM[i - 1][j] ) {
					SM[i][j] = SM[i - 1][j] ;
					tmp = 1;
				}
			}
			if (SM[i][j] < SM[i - 1][j]+gap) {
				SM[i][j] = SM[i - 1][j]+gap;
				tmp = 1;
			}
			int tmp1 = SM[i][j];
			if (a != -1 && b != -1) {
				SM[i][j] = max(SM[i][j], SM[i - 1][j - 1]+matrix[a][b]);
			}
			if ((a == -1 && b != -1)||(i!=-1&&j==-1)) {
				SM[i][j] = max(SM[i][j], SM[i - 1][j - 1] + gap);
			}
			if (a == -1 && b == -1) {
				SM[i][j] = max(SM[i][j], SM[i - 1][j - 1]);
			}
			if (SM[i][j] != tmp1) {
				tmp = 2;
			}
			if (tmp == 0) {
				ins1[i][j] = ins1[i][j - 1] + 1;
				ins2[i][j] = ins2[i][j - 1] + 1;
				S1[i][j] = S1[i][j - 1];
				S2[i][j] = S2[i][j - 1];
				if (ins1[i][j - 1] > S1[i][j-1].size()) {
					S1[i][j].insert(S1[i][j].size(), "-");
				}
				else {
					S1[i][j].insert(ins1[i][j - 1], "-");
				}
			}
			if (tmp == 1) {
				ins1[i][j] = ins1[i - 1][j] + 1;
				ins2[i][j] = ins2[i - 1][j] + 1;
				S2[i][j] = S2[i - 1][j];
				S1[i][j] = S1[i - 1][j]; 
				if (ins2[i - 1][j] > S2[i - 1][j].size()) {
					S2[i][j].insert(S2[i - 1][j].size(), "-");
				}
				else {
					S2[i][j].insert(ins2[i - 1][j], "-");
				}
			}
			if (tmp == 2) {
				S1[i][j] = S1[i - 1][j - 1];
				S2[i][j] = S2[i - 1][j - 1];
				ins1[i][j] = ins1[i - 1][j-1] + 1;
				ins2[i][j] = ins2[i - 1][j-1] + 1;
			}
		}
	}

	return make_tuple(SM[s1.size()][s2.size()],S1[s1.size()][s2.size()],S2[s1.size()][s2.size()]);
}
int main() {
	string copia=  "ILDFHEKLLHPGIQKTTKLFGETYYFPNSQLLIQNIINECSICNLAK";
	string MMULV = "LLDFLLHQLTHLSFSKMKALLERSHSPYYMLNRDRTLKNITETCKACAQVN";
	string HTLV =  "LQLSPAELHSFTHCGQTALTLQGATTTEASNILRSCHACRGGN";
	string RSV =   "YPLREAKDLHTALHIGPRALSKACNISMQQAREVVQTCPHCNSA";
	string MMTV =  "IHEATQAHTLHHLNAHTLRLLYKITREQARDIVKACKQCVVAT";
	string SMRV =  "LESAQESHALHHQNAAALRFQFHITREQAREIVKLCPNCPDWGS";
	string list[6] = { copia,MMULV,HTLV,RSV,MMTV,SMRV };
	int gap = -2;
	int max = INT_MIN;
	int argmax = 0;
	int Sc;
	tuple<int, string, string> test;
	//(1)全てのシーケンスの類似度を測定して類似度の和を最大化するシーケンスを選ぶ
	for (int i = 0; i < 6; i++) {
		Sc = 0;
		for (int j = 0; j < 6; j++) {
			if (i == j)continue;
			//cout << list[i] << " "<<list[j] << endl;
			test=(score(list[i], list[j], gap));
			Sc += get<0>(test);
			//cout << get<0>(test) << endl;
			//cout << get<1>(test) << endl;
			//cout << get<2>(test) << endl;
		}
		if (max < Sc) {
			max = Sc;
			argmax = i;
			//cout << max << endl;
		}
		
	}
	//(2) (1)で選んだシーケンスと他のシーケンスのアライメントを行う。
	vector<tuple<int, string, string>> list1(0);
	for (int i = 0; i < 6; i++) {
		if (i == argmax)continue;
		test = score(list[argmax], list[i],gap);
		//cout << get<0>(test) << endl;
		//cout << get<1>(test) << endl;
		//cout << get<2>(test) << endl;
		get<0>(test)= 0;
		
		list1.push_back(test);
	}
	//(3) (1)で選んだシーケンスを基にして(2)のアライメントからmultiple alignment を行う
	string Slist[6];
	for (int i = 0; i < 5;i++ ) {
		tuple<int, string, string > Tup= list1[i];//loopを確認するため
		while (get<0>(list1[i]) < get<1>(list1[i]).size()) {
			bool tmp = true;
			char ins = 0;
			for (int j = 0; j < 5; j++) {
				int index = get<0>(list1[j]);
				ins = get<1>(list1[j])[index];

				if (ins == '-') {
					tmp = false;
				}
				
			}
			if (tmp) {
				Slist[argmax].insert(Slist[argmax].size(), { ins });
				for (int j = 0; j < 5; j++) {
					int index =  get<0>(list1[j]);
					ins = get<2>(list1[j])[index];
					int index1 = j;
					if (j >= argmax)index1 = index1 + 1;
					if (index >= get<2>(list1[j]).size())ins = '-';
					Slist[index1].insert(Slist[index1].size(), { ins });
					get<0>(list1[j])= get<0>(list1[j])+1;
				}
				get<0>(Tup) = get<0>(Tup) + 1;
			}
			else {
				Slist[argmax].insert(Slist[argmax].size(), "-");
				for (int j = 0; j < 5; j++) {
					int index = get<0>(list1[j]);
					ins = get<1>(list1[j])[index];
					int index1 = j;
					if (j >= argmax)index1 = index1 + 1;
					if (index >= get<1>(list1[j]).size())ins = '-';
					if (ins == '-') {
						Slist[index1].insert(Slist[index1].size(), {get<2>(list1[j])[index]});
						get<0>(list1[j])= get<0>(list1[j]) +1;
					}
					else {
						Slist[index1].insert(Slist[index1].size(), "-");

					}				
				}
			}
			
		}
	}
	for (int i = 0; i < 6; i++) {
		cout << Slist[i] << endl;
	}
}
