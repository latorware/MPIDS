#include <iostream>
#include <vector>
#include <set>
#include <math.h>

using namespace std;

void leer_grafo(vector<set<int>>& grafo) {
	
	
	int u, v;
    while (cin >> u >> v) {
        grafo[u - 1].insert(v - 1);
        grafo[v - 1].insert(u - 1);
    }

}

void escribir_grafo(const vector<set<int>>& grafo) {
	for (int i = 0; i < grafo.size(); ++i) {
	   cout << i << ":";
	   for (int j : grafo[i]) {
	   		cout << " " << j; 
		}
	   cout << endl;
	}
}

bool is_minimal(const vector<set<int>>& grafo, const set<int>& set) {
	bool is = true;	
	
	for (int i : set) if(is)	{
	bool out = true;
	   for (int j : grafo[i]) if (out) {
	   int size = grafo[j].size();
	   size = ceil(size/2.0);
	   
	   for (int k : grafo[j]) {
		  if (set.find(k) != set.end() and k != i) --size; 
		  
	   }
	   if (size > 0) out = false;
	}	
	if (out) is = false;
	}
	
	return is;
	
}

bool is_dominant(const vector<set<int>>& grafo, const set<int>& set) {
	bool is = true;
	for (int i = 0; i < grafo.size() && is; ++i) {
	   int size = grafo[i].size();
	   size = ceil((size/2.0));
	   
	   for (int j : grafo[i]) {
		  bool found = false;
		  for (int k : set) if(!found){
			 if (j == k) {
			   --size;
			   found = true;
				 
			 }
		  }
		   
	   }
	   if (size > 0) {
		   is = false;  
	   
	   }
	}
	return is;
}

int main() {
	int n, m;
	cin >> n >> m;
	
	vector<set<int>> grafo(n);
	leer_grafo(grafo);
	escribir_grafo(grafo);
	set<int> set = {0, 2, 4, 5, 6, 9};
	
	
	if (is_dominant(grafo, set) ) {
	  cout << "Dominant: True" << endl;
	  if (is_minimal(grafo, set)) cout << "Minimal: True" << endl;
	  else cout << "Minimal: False" << endl;
	} else cout << "Dominant: False" << endl;
	
	
	
// 	cout << " I am ceiling: "  << ceil((3/2.0)); 
	
}
