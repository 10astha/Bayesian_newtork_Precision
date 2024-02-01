#include<bits/stdc++.h>

using namespace std;
using ll = long long;
using ull = unsigned long long;

#define dbg(x) cout << x << "\n";
#define pr(x) cout << fixed << setprecision(16) << x << "\n";

// Taken from Um_NIK 
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

ll myRand(ll B) {
	return (ull)rng() % B;
}
clock_t startTime;
double getCurrentTime() {
	return (double)(clock() - startTime) / CLOCKS_PER_SEC;
}



class Graph_Node{
public:
//************************************************************************************************************************//
	string Node_Name;
    unordered_map<int,int> Children;
    unordered_map<string,int> rnkkk;
	vector<string> Parents;
	int nvalues;
	vector<string> values;
	vector<double> CPT;
	vector<double> Ct;
//************************************************************************************************************************//
    Graph_Node(string name,int n,vector<string> vals) { Node_Name = name; nvalues = n; values = vals;
        for(int i = 0; i < vals.size(); i++) rnkkk[vals[i]] = i;
     }
//************************************************************************************************************************//
	string get_name() { return Node_Name;}
//************************************************************************************************************************//
	vector<int> get_children() { 
        vector<int> vv; 
        for(auto it = Children.begin(); it != Children.end(); it++) {
            vv.push_back(it -> first);
        }   
        return vv; 
    }
//************************************************************************************************************************//
	vector<string> get_Parents() { return Parents; }
//************************************************************************************************************************//
	vector<double> get_CPT() { return CPT; }
//************************************************************************************************************************// 
	int get_nvalues() { return nvalues; }
//************************************************************************************************************************//
	vector<string> get_values() { return values; }
//************************************************************************************************************************//
	void set_CPT(vector<double>& new_CPT) { CPT.clear(); CPT=new_CPT; }
//************************************************************************************************************************//
	void reset_Ct(){
		Ct.resize(CPT.size());
		for(int i = 0; i < CPT.size(); i++) Ct[i] = 0.0035;
	}
//************************************************************************************************************************//
	void set_Parents(vector<string> Parent_Nodes) { Parents.clear(); Parents = Parent_Nodes; }
//************************************************************************************************************************//
	void add_child(int new_child_index) { 
        Children[new_child_index] = 1;
	}
//************************************************************************************************************************//
	int get_rank(string given_value){
        return rnkkk[given_value];
	}
//************************************************************************************************************************//
	float get_probab(int i){ return CPT[i]; }
//************************************************************************************************************************//	
	void init_pct(){
		double p = (1.0)/(1.0 * nvalues);
		for(int i = 0 ; i < CPT.size(); i++) {
			CPT[i] = p;
		}
		return;
	}
//************************************************************************************************************************//
};

class network{
public:
//************************************************************************************************************************//
	list <Graph_Node> Pres_Graph;
	vector<string> var_names;
	vector<vector<string>> missing_data;
    vector<vector<vector<int>>> arr;
	int addNode(Graph_Node node) { Pres_Graph.push_back(node); return 0; }
//************************************************************************************************************************//
	int netSize() { return Pres_Graph.size();  }
//************************************************************************************************************************//
	int get_index(string var_name){ list<Graph_Node>::iterator listIt; int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++) { if(listIt->get_name() == var_name) return count; count++; }
        return -1; 
	}
//************************************************************************************************************************//
    list<Graph_Node>::iterator get_nth_node(int n) {
       list<Graph_Node>::iterator listIt; int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++) { 
			if(count==n) return listIt; count++; } 
			return listIt;  
	}
//************************************************************************************************************************//
    list<Graph_Node>::iterator search_node(string val_name) {
        list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++) {
            if(listIt->get_name().compare(val_name)==0) return listIt;
        }
        // cout<<"node not found\n";
        return listIt;
    }
//************************************************************************************************************************//
//  laplace add one smoothing.
	void reset_Ct_net(){
		auto itr = Pres_Graph.begin();
		for( ; itr != Pres_Graph.end(); itr++) {
			for(int i = 0; i < itr->Ct.size(); i++) {
				itr -> Ct[i] = 0.0035;
			}
		}
	}
//************************************************************************************************************************//
	int get_cpt_index(list<Graph_Node>::iterator &child, vector<string> &values){
		vector<string> par = child -> get_Parents();
		int n = par.size();
		int curr_prod = 1;
		int ans = 0;
		for(int i = n - 1; i >= 0; i--){
			list<Graph_Node>::iterator itr = search_node(par[i]);  // current parent.
			int idx = get_index(par[i]);   // it's index.
			int rnk = itr->get_rank(values[idx]); // it's value's rank.
			int sz = itr->get_nvalues(); // it's number of values.
			ans += rnk * curr_prod; // 
			curr_prod *= sz; //
		}
		int idx = get_index(child -> get_name());  // from name to index.
		int rnk = child -> get_rank(values[idx]);  // need to verify this indeed gives correct value.
		ans += rnk * curr_prod;
		return ans;
	}
//************************************************************************************************************************//
	double cnd_probab(list<Graph_Node>::iterator &child, vector<string> &data) {
		int index = get_cpt_index(child, data);
		double p = child -> get_probab(index);
		return p;
	}
//************************************************************************************************************************//
	void init() {
		int n = netSize();
        int cnt = 0;
		for(auto itr = Pres_Graph.begin(); itr != Pres_Graph.end(); itr++) {
			int x = itr -> nvalues;
			itr->reset_Ct();
			itr->init_pct();
			cnt++;
		}
        arr.resize(missing_data.size());
        for(int i = 0; i < missing_data.size(); i++){
            int k = stoi(missing_data[i].back());
			if(k == -1){ // complete data.
				arr[i].resize(1);
				arr[i][0].resize(n);
			}
			else {
				auto itr = get_nth_node(k);
            	k = itr->nvalues;
            	arr[i].resize(k);
            	for(int j = 0; j < k; j++) {
                	arr[i][j].resize(n);
            	}
			}
        }
		return;
	}
//************************************************************************************************************************//
	double PmB(list<Graph_Node>::iterator &itr, int itr_idx, vector<string> &data){
		string val = data[itr_idx];
		double result , answer;  
		vector<int> child = itr -> get_children();  
		int total_values = itr -> nvalues;
		vector<double> Probab(total_values);
		for(int i = 0; i < total_values; i++) {
			data[itr_idx] = itr->values[i];
            answer = 0.0;
			for(auto &v : child) {
				auto utr = get_nth_node(v);
				result = cnd_probab(utr, data);
                answer += log(result);
			}
			result = cnd_probab(itr, data);
            answer += log(result);
            Probab[i] = exp(answer);
		}
		answer = 0.0;
		for(int i = 0; i < total_values; i++){
			answer += Probab[i];
		}
        answer = exp(log(Probab[itr->get_rank(val)]) - log(answer));
		data[itr_idx] = val;
        return answer;
	}
//************************************************************************************************************************//
	void add_data(int i, int val_rnk, float weight) {
		if(val_rnk == -1) {
			int cnt = 0;
			for(auto itr = Pres_Graph.begin(); itr != Pres_Graph.end(); itr++) {
				int idx = arr[i][0][cnt];
				itr ->Ct[idx] += weight;
				cnt++;
			}
			return;
		}
        int cnt = 0;
        for(auto itr = Pres_Graph.begin(); itr != Pres_Graph.end(); itr++) {
            int idx = arr[i][val_rnk][cnt];
			itr->Ct[idx] += weight;
            cnt++;
        }
		return;
	}
//************************************************************************************************************************//
	void add_data_1(int i ,int rnk, vector<string> data) {
		if(rnk == -1) { // it is a complete data.
			int cnt  = 0;
			for(auto itr = Pres_Graph.begin(); itr != Pres_Graph.end(); itr++) {
				arr[i][0][cnt] = get_cpt_index(itr, data);
				cnt++;
			}
			return;
		}
        int cnt = 0;
        for(auto itr = Pres_Graph.begin(); itr != Pres_Graph.end(); itr++) {
			arr[i][rnk][cnt] = get_cpt_index(itr, data);
            cnt++;
        }
		return;
	}
//************************************************************************************************************************//
    void pre_process(vector<vector<string>> &Data){
        vector<string> d;
        for(int i = 0; i < Data.size(); i++){
            d = Data[i];
			int j = stoi(d.back());
			if(j == -1){
				add_data_1(i, -1, d);
			}		
			else{
                auto itr = get_nth_node(j);
				vector<string> vb = itr->values;
                int cntrrr = 0;
				for(auto v : vb){
					d[j] = v;
					float w = PmB(itr, j, d);
                    add_data_1(i,cntrrr, d);
                    cntrrr++;
				}
			}
		}
    }
//************************************************************************************************************************//	
	void func(vector<vector<string>> &Data){
		vector<string> d;
		for(int i = 0; i < Data.size(); i++){
			d = Data[i];
			int j = stoi(d.back());
			if(j == -1){
				add_data(i, -1, 1.0);
			}
			else{
				auto itr = get_nth_node(j);
				vector<string> vb = itr->values;
                int cntk = 0;
				for(auto v : vb){
					d[j] = v;
					float w = PmB(itr, j, d);
					add_data(i,cntk, w);
                    cntk++;
				}
			}
		}
	}
//************************************************************************************************************************//s
	void pct_from_ct(int idx_v) {
		auto itr = get_nth_node(idx_v);
		vector<string> par = itr -> get_Parents();
		int P = 1;
		int res;
		int s0 = itr -> get_nvalues();
		for(auto &x : par){
			auto utr = search_node(x);
			res = utr -> get_nvalues();
			P *= res;
		}
		for(int i = 0; i < P; i++) {
			double sum = 0.0;
            int t = i;
			for(int rk = 0; rk < s0; rk++) {
				sum += itr->Ct[t];
                t += P;
			}
            t = i;
			for(int rk = 0; rk < s0; rk++) {
				itr->CPT[t] = (itr->Ct[t])/(sum);
                t += P;
			}
		}
	}
//************************************************************************************************************************//
	void update_pct(){
		for(int i = 0; i < netSize(); i++) {
			pct_from_ct(i);
		}
	}	
//************************************************************************************************************************//
	void print(){
		auto itr = Pres_Graph.begin();
		itr++;
		for(int i = 0; i < itr->CPT.size(); i++){
			pr(itr->CPT[i]);
		}
	}
};

double TT = 113.0;

network read_network(string f1, string f2, string f3) {
	network Alarm;
	string line;
	int find = 0;
  	ifstream myfile1(f1); 
  	string temp;
  	string name;
  	vector<string> values;
	if (myfile1.is_open()) {
		while(!myfile1.eof()) {	
			stringstream ss;
      		getline (myfile1,line);
      		ss.str(line);
     		ss >> temp;
			if(temp == "variable") {
     				ss >> name;
     				getline (myfile1,line);
     				stringstream ss2;
     				ss2.str(line);
     				for(int i = 0; i < 4; i++) {ss2 >> temp;}
     				values.clear();
     				while(temp != "};") { values.push_back(temp); ss2 >> temp;}
     				Graph_Node new_node(name, values.size(), values);
					Alarm.addNode(new_node);
     		}
			else if(temp == "probability") {        
					ss >> temp;
     				ss >> temp;
                    list<Graph_Node>::iterator listIt;
                    list<Graph_Node>::iterator listIt1;
     				listIt = Alarm.search_node(temp);
                    int index = Alarm.get_index(temp);
                    ss >> temp;
                    values.clear();
					while(temp != ")") {
                        listIt1 = Alarm.search_node(temp);
                        listIt1 -> add_child(index);
     					values.push_back(temp);
     					ss >> temp;
    				}
                    listIt -> set_Parents(values);
    				getline (myfile1,line);
     				stringstream ss2;
     				ss2.str(line);
     				ss2 >> temp;
     				ss2 >> temp;
     				vector<double> curr_CPT;
                    string::size_type sz;
     				while(temp != ";") {
     					curr_CPT.push_back(atof(temp.c_str()));	
     					ss2 >> temp;
    				}
                    listIt->set_CPT(curr_CPT);
			}
            else{}
    	}
    	if(find==1)
    	myfile1.close();
  	}

	ifstream myfile2(f2);
	int n = Alarm.netSize();
	while(!myfile2.eof()){
		    string line;
			stringstream ss;
      		getline (myfile2,line);
      		ss.str(line);
			vector<string> lien;
			int n = Alarm.netSize();
			string tmp;
			int cntr = -1;
			for(int i = 0; i < n; i++) {
				ss >> tmp;
				if(tmp == "\"?\"") {cntr = i;}
				lien.push_back(tmp);
			}
			string cnt = to_string(cntr);
			lien.push_back(cnt);
			(Alarm.missing_data).push_back(lien);
	}


	Alarm.init();
    Alarm.pre_process(Alarm.missing_data);


	double curr_time = getCurrentTime();
	while(curr_time < TT) {
		Alarm.func(Alarm.missing_data);
		Alarm.update_pct();
		Alarm.reset_Ct_net();
		curr_time = getCurrentTime();
	}


    ifstream myfile3(f1);
    string temp_input;
    ofstream out;
    string word;
    out.open(f3);
    if (!myfile3.is_open()){ return Alarm;}
    else{
        while (!myfile3.eof()){
            stringstream ss;
            getline (myfile3,line);
            ss.str(line);
            ss >> word;
            if(word=="probability"){
                ss >> word;
                ss >> word;
                int word_index = Alarm.get_index(word);
                out << line << endl;
                getline(myfile3,line);
                out<< "    table ";  
                auto itr = Alarm.get_nth_node(word_index);
                int cptsize = itr->CPT.size();
                vector<double> cpt = itr->get_CPT();
                for(int i = 0; i < cptsize; i++){
                    out << fixed << setprecision(4) << cpt[i] << " ";
                }
                out << ";" << endl;
            }
            else if(line != ""){ out << line << endl;}
            else{out << line;}
        }
    }
    return Alarm;
}

int main(int argc, char* argv[]) {
	string x1 = argv[1];
	string x2 = argv[2];
    string x3 = "solved_alarm.bif";
	network Alarm;
	Alarm = read_network(x1, x2, x3);
	int n = Alarm.netSize();
	auto end = getCurrentTime();
	pr(end);
}