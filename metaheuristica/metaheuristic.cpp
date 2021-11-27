/***************************************************************************
    metaheuristic.cpp 
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

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// computing time limit for each application of the metaheuristic
double time_limit = 3200.0;

// number of applications of the metaheuristic
int n_apps = 1;
int pob = 10;
int gen = 10;

// dummy parameters as examples for creating command line parameters 
// (see function read_parameters(...))
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
        if (strcmp(argv[iarg],"-i") == 0) inputFile = argv[++iarg];
        // reading the computation time limit 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-t") == 0) time_limit = atoi(argv[++iarg]); 
        // reading the number of applications of the metaheuristic 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps") == 0) n_apps = atoi(argv[++iarg]); 
        // example for creating a command line parameter 
        // param1 -> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1") == 0) {
            dummy_integer_parameter = atoi(argv[++iarg]); 
        }
        // example for creating a command line parameter 
        // param2 -> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2") == 0) {
            dummy_double_parameter = atof(argv[++iarg]);
        }
        else if (strcmp(argv[iarg],"-pob") == 0) {
            pob = atoi(argv[++iarg]); 
        }
        else if (strcmp(argv[iarg],"-gen") == 0) {
            gen = atoi(argv[++iarg]); 
        }
        iarg++;
    }
}

//Dona un int del 0 al x-1
int random(int x) {
    float n = rnd->next();
    return (int)floor(n*x) % x;
}




bool is_dominant(const vector<set<int>>& grafo, const set<int>& set) {
	bool is = true;
	for (int i = 0; i < grafo.size() && is; ++i) {
	   int size = grafo[i].size();
	   size = ceil((size/2.0));
	   
	   for (int j : grafo[i]) {
		  bool found = false;
		  for (int k : set) {
			if (!found and j == k) {
			   --size;
			   found = true;
			}
		  }
	   }
	   if (size > 0) is = false;
	}
	return is;
}

set<int> dominantSet(vector<bool> scenario) {
    set<int> result;
    for (int i = 0; i < scenario.size(); ++i) {
        if (scenario[i]) result.insert(i);
    }
    return result;
}


vector<bool> randomSolution() {
    vector<bool> sol(n_of_nodes);
    set<int> dominant = set<int>();
    for (int i = 0; i < n_of_nodes; ++i) {
        if (rnd->next() <= 0.9) {
            sol[i] = true; 
            dominant.insert(i);
        }
        else sol[i] = false;
    }
    if (is_dominant(neighbors, dominant)) return sol;
    return randomSolution();
}

float fitness(vector<bool> sol) {
    float sub = 0;
    for (int i = 0; i < sol.size(); ++i) {
    	if (sol[i]) ++sub;
    }
    return (float)n_of_nodes - sub;
}


//el bool t es un parametre per referencia que ha de posar true si sha trobat dues combinacions valides entre dos solucions
pair<vector<bool>,vector<bool>> randCombination(const vector<bool>& v1, const vector<bool>& v2) {
    vector<bool> v1_aux = v1;
    vector<bool> v2_aux = v2;
    bool d1 = false;
    bool d2 = false;
    int j = 0, r;
    cout << "recombination" << endl;
    while ((j < n_of_nodes*0.75) and (not d1 or not d2)) {
        r = random(v1.size()-1);
        if (d1) {
            v2_aux = v2;
            for (int i = 0; i <= r; ++i) {
                v2_aux[i] = v1[i];
            }
            d2 = is_dominant(neighbors, dominantSet(v2_aux));
        }
        else if (d2) {
            v1_aux = v1;
            for (int i = 0; i <= r; ++i) {
                v1_aux[i] = v2[i];
            }
            d1 = is_dominant(neighbors, dominantSet(v1_aux));
        }
        else {
            v1_aux = v1;
            v2_aux = v2;
            for (int i = 0; i <= r; ++i) {
                v1_aux[i] = v2[i];
                v2_aux[i] = v1[i];
            }
            d1 = is_dominant(neighbors, dominantSet(v1_aux));
            d2 = is_dominant(neighbors, dominantSet(v2_aux));
        }
        ++j;
    }
    //cout << "recombination2" << endl;
    if (d1 and d2) return make_pair(v1_aux,v2_aux);
    if (d1) return make_pair(v1_aux,v2);
    if (d2) return make_pair(v1,v2_aux);
    return make_pair(v1, v2);
}

vector<bool> randMutation(const vector<bool>& v) {
    vector<bool> v_aux = v; 
    set<int> s = dominantSet(v);
    int j = 0;
    bool d = false;
    while (j < n_of_nodes*0.75 and not d) {
        int r = random(v.size());
        v_aux[r] = not v_aux[r];
        if (v_aux[r]) s.insert(r);
        else s.erase(r);
        ++j;
        d = is_dominant(neighbors, s);
    }
    //cout << "mutation" << endl;
    if (d) return v_aux;
    return v;
}



void addInOrderLim(vector<pair<vector<bool>,float>> &best, vector<bool> next) {
	float ftn = fitness(next);
	int pos = -1;
	for (int i = 0; i < best.size() and pos == -1; ++i) {
		if (ftn > best[i].second) pos = i;
	}
	
	for (int j = best.size()-2; j >= pos and pos != -1; --j) {
		best[j+1] = best[j];
	}
	if (pos != -1) {
		pair<vector<bool>,float> aux(next,ftn);
		best[pos] = aux;		
	}
}

void addInOrder(vector<pair<vector<bool>,float>> &best, vector<bool> next) {
	float ftn = fitness(next);
	int pos = -1;
	for (int i = 0; i < best.size() and pos == -1; ++i) {
		if (ftn > best[i].second) pos = i;
	}
	for (int j = best.size()-2; j >= pos and pos != -1; --j) best[j+1] = best[j];

    pair<vector<bool>,float> aux(next,ftn);
	if (pos != -1) {best[pos] = aux;}
    else best.push_back(aux);
}


/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. A random number 
    // between 0 and 1 is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // vectors for storing the result and the computation time 
    // obtained by the <n_apps> applications of the metaheuristic
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if (not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< set<int> >(n_of_nodes);
    int u, v;
    while (indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    // main loop over all applications of the metaheuristic
    for (int na = 0; na < n_apps; ++na) {

        // the computation time starts now
        Timer timer;

        // Example for requesting the elapsed computation time at any moment: 
        // double ct = timer.elapsed_time(Timer::VIRTUAL);

        cout << "start application " << na + 1 << endl;

        // HERE GOES YOUR METAHEURISTIC
        vector<pair<vector<bool>,float>> population = vector<pair<vector<bool>,float>>(pob);
        
        
        for (int i = 0; i < pob and timer.elapsed_time(Timer::VIRTUAL) < time_limit; ++i) {
            cout << "randomSol" << endl;
            addInOrder(population, randomSolution());
            cout << "out of randomSol" << endl;
        }
        for (int g = 0; g < gen and timer.elapsed_time(Timer::VIRTUAL) < time_limit; ++g) {
            //Conservem un 20% de les millors solucions
            vector<vector<bool>> new_population = vector<vector<bool>>(pob);
            for (int cross = 0; cross < floor(pob*0.4) and timer.elapsed_time(Timer::VIRTUAL) < time_limit; ++cross) {
            	pair<vector<bool>,vector<bool>> aux;
            	int a, b;
                a = random(pob);
        	    while (rnd->next() < 0.4 and timer.elapsed_time(Timer::VIRTUAL) < time_limit) a = random(pob);
        	    b = random(pob);
        	    while (rnd->next() < 0.4 and a != b and timer.elapsed_time(Timer::VIRTUAL) < time_limit) b = random(pob);
        	    if (timer.elapsed_time(Timer::VIRTUAL) < time_limit)  {
                    aux = randCombination(population[a].first,population[b].first);
                    new_population[2*cross] = aux.first;
                    new_population[2*cross+1] = aux.second;
                }
            }

            for (int k = 0; k < floor(pob*0.4) and timer.elapsed_time(Timer::VIRTUAL) < time_limit; ++k) {
                if (rnd->next() < 0.1) new_population[k] = randMutation(new_population[k]);
            }

            for (int k = 0; k < floor(pob*0.4) and timer.elapsed_time(Timer::VIRTUAL) < time_limit; ++k) {
                addInOrderLim(population, new_population[k]);
            }

        }
        double ct = timer.elapsed_time(Timer::VIRTUAL);

        
        
        // Whenever the best found solution is improved, first take the 
        // computation time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: 
        cout << "value " << n_of_nodes - population[0].second;
        cout << "\ttime " << ct << endl;
        //for (int i = 0; i < population[0].first.size(); ++i) if (population[0].first[i]) cout <<"-" <<i <<"-";
        //

        // Store the value of the new best found solution in vector results: 
        results[na] = n_of_nodes - population[0].second;
        //
        // And store the current computation time (that is, the time measured 
        // at that moment and stored in variable "ct") in vector times: 
        times[na] = ct;
        //
        // Stop the execution of the metaheuristic 
        // once the time limit "time_limit" is reached.

        cout << "end application " << na + 1 << endl;
    }

    // calculating the average of the results and computation times, 
    // and their standard deviations, and write them to the screen
    double r_mean = 0.0;
    int r_best = std::numeric_limits<int>::max();
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        if (int(results[i]) < r_best) r_best = int(results[i]);
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    double rsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd/double(results.size());
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    tsd = tsd/double(results.size());
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }
    // printing statistical information
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

