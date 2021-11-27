/***************************************************************************
    greedy.cpp 
    (C) 2021 by C. Blum & M. Blesa
    
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
#include "Random.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>

#include <algorithm>

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// dummy parameters as examples for creating command line parameters 
// see function read_parameters(...)
int dummy_integer_parameter = 0;
int dummy_double_parameter = 0.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        
        // example for creating a command line parameter param1 
        //-> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) 
            dummy_integer_parameter = atoi(argv[++iarg]); 
            
        // example for creating a command line parameter param2 
        //-> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) 
            dummy_double_parameter = atof(argv[++iarg]);  
            
        iarg++;
    }
}


bool ordena_es_major(int f, int g) { //per ordenar veins de major a menor
    //cout << f << ' ' << g << endl; 
    if (neighbors[f].size() >= neighbors[g].size()) return true;
    else return false; 
}


bool es_dominant(set<int>& c_dominant, vector< set<int> >& neighbors) { //NOMES DOMINANT (no influencia positiva)
    for (int i = 0; i < neighbors.size(); i++) {
        if (    !(c_dominant.find(i) != c_dominant.end())     ) {
            bool trobat_almemys_un = false; 
            set<int>::iterator it; 
            for (it = neighbors[i].begin(); (it != neighbors[i].end() && (!trobat_almemys_un)); ++it) {
                if ( (c_dominant.find(*it) != c_dominant.end())     ) trobat_almemys_un = true; 
            }
            if (!trobat_almemys_un) return false; 
        }
    }
    return true; 
}


bool inf_positiva(set<int>& c_inf_positiva, vector<set<int> >& neighbors) {
    for (int i = 0; i < neighbors.size(); i++) {
        int meitat = ceil(neighbors[i].size()/2.0); 
        bool almenys_meitat = false; 
        int contador = 0; 
        set<int>::iterator it; 
        for (it = neighbors[i].begin(); (it != neighbors[i].end()) && (!almenys_meitat); it++) {
            if ( (c_inf_positiva.find(*it) != c_inf_positiva.end())) contador++; 
            if (contador >= meitat) almenys_meitat = true; 
        }
        if (!almenys_meitat) return false; 

    }
    return true; 
}
/**
void escribir_grafo( vector<set<int> >& grafo) {
	for (int i = 0; i < grafo.size(); ++i) {
	   cout << i << ":";
       if (grafo[i].empty()) cout << "TRUE" << endl; 
	   for (int j : grafo[i]) {
	   		cout << " " << j; 
		}
	   cout << endl;
	}
}
*/

void MergeSortedIntervals(vector<int>& v, int s, int m, int e) {
	
    // temp is used to temporary store the vector obtained by merging
    // elements from [s to m] and [m+1 to e] in v
	vector<int> temp;

	int i, j;
	i = s;
	j = m + 1;

	while (i <= m && j <= e) {

		if (ordena_es_major(v[i], v[j])) {
			temp.push_back(v[i]);
			++i;
		}
		else {
			temp.push_back(v[j]);
			++j;
		}

	}

	while (i <= m) {
		temp.push_back(v[i]);
		++i;
	}

	while (j <= e) {
		temp.push_back(v[j]);
		++j;
	}

	for (int i = s; i <= e; ++i)
		v[i] = temp[i - s];

}

// the MergeSort function
// Sorts the array in the range [s to e] in v using
// merge sort algorithm
void MergeSort(vector<int>& v, int s, int e) {
	if (s < e) {
		int m = (s + e) / 2;
		MergeSort(v, s, m);
		MergeSort(v, m + 1, e);
		MergeSortedIntervals(v, s, m, e);
	}
}



int hs (int a, const set<int>& result) {
    int meitatdeg = ceil(neighbors[a].size()/2.0); 
    int veins_dominant = 0; 
    for (int i : neighbors[a]) {
        if ((result.find(i) != result.end())) veins_dominant++; 
    }
    return (meitatdeg - veins_dominant); 
}

int coverdegree (int a, const set<int>& result) {
    int contador = 0; 
    for (int i : neighbors[a]) {
        if (hs(i, result) > 0) contador++; 
    }
    return contador; 
}


/************
Main function
*************/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. 
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // variables for storing the result and the computation time 
    // obtained by the greedy heuristic
    double results = std::numeric_limits<int>::max();
    double time = 0.0;

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< set<int> >(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    // the computation time starts now
    Timer timer;

    // Example for requesting the elapsed computation time at any moment: 
    // double ct = timer.elapsed_time(Timer::VIRTUAL);

    // HERE GOES YOUR GREEDY HEURISTIC
    // When finished with generating a solution, first take the computation 
    // time as explained above. Say you store it in variable ct.
    // Then write the following to the screen: 
    // cout << "value " << <value of your solution> << "\ttime " << ct << endl;

    //escribir_grafo(neighbors); 
    set<int> resultat = set<int> (); //Conjunt dominant dinfluencia positiva
    //set<int> inverseresult = set<int> (); 


    //ORDENAR PER NOMBRE VEINS
    vector<int> veins = vector<int> (); //veins[i] sera el vertex ordenat per mes a menys arestes

    for (int i = 0; i < neighbors.size(); i++) {
        veins.push_back(i); 
        //inverseresult.insert(i); 
    }
    /**
    cout << "IMPRIMINT RESULTAT: "; 
    for (int i : resultat) cout << i << ' '; 
    cout << endl; 
    cout << "IMPRIMINT INVERSE: "; 
    for (int i : inverseresult) cout << i << ' '; 
    cout << endl; 
    //cout << "VEINS SIZE " << veins.size() << endl; 
    //sort (veins.begin(), veins.end(), ordena_es_major); 
    */

    MergeSort(veins, 0, veins.size()-1); 

    cout << "time after merge sort " << timer.elapsed_time(Timer::VIRTUAL) << endl; 
    //for (int i : veins) cout << i << ' '; 
    //cout << endl; 
    //JA TENIM VEINS (veins[i]) ORDENATS DE MAJOR A MENOR PER ARESTES

    bool dominant_positiu = false; //diu si es dominant (nomes dominant, no dominant de influencia positiva)
    int d_a_partir_de = 0; 
    //fem fins que el tinguem dominant
    for (int i = 0; (i < veins.size()); i++) {
        int hsi = hs(veins[i], resultat);
        if ( hsi > 0 ) {
            //cout << "ENTRA EN HS > 0 " << veins[i]; 
            //cout << "ENTRANDO EN HSI > 0     " << veins[i] << endl; 
            for (int j = 0; j < hsi; j++) {
                //cout << "ENTRANDO EN BUCLE 0 HSI     "<< j << endl; 
                bool es_primer = true; 
                bool almenys_un = false; 
                pair <int, int> maximo; //first = vetrtice_maximo   second = cover degree
                for (int k : neighbors[veins[i]]) {
                    //cout << "ENTRANDO EN BUCLE k NEIGHOURS     "<< k <<endl; 
                    //cout << k << ' '; 
                    if (resultat.find(k) == resultat.end()) {
                        if (!almenys_un) almenys_un = true; 
                        if (es_primer) {
                            maximo.first = k; 
                            maximo.second = coverdegree(k, resultat); 
                            es_primer = false; 
                        }
                        else {
                            int actual = coverdegree(k, resultat); 
                            if (actual >= maximo.second) {
                                maximo.first = k; 
                                maximo.second = actual; 
                            }

                        }
                    }
                }
                //cout << endl; 
                //cout << endl; 
                if (almenys_un) {
                    resultat.insert(maximo.first); 
                }

            }
            //cout << endl; 
        }

        //cout << endl;
    }

    /**
    cout << "IMPRIMINT RESULTAT: "; 
    for (int i : resultat) cout << i << ' '; 
    cout << endl; 
    cout << "IMPRIMINT INVERSE: "; 
    for (int i : inverseresult) cout << i << ' '; 
    cout << endl; 
    */

    if (true) {
        double ct = timer.elapsed_time(Timer::VIRTUAL);
        cout << "value " << resultat.size() << "\ttime " << ct << endl;
        //for(int i : resultat) cout<< i << endl;
    }
    else {
        cout << "ERROR. " << endl; 
    }


}

