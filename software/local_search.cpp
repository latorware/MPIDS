/***************************************************************************
    local_search.cpp
    (C) 2021 by C.Blum & M.Blesa
    
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
#include <math.h>

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;
vector<int> veins_dominants; 

// string for keeping the name of the input file
string inputFile;

// number of applications of local search
int n_apps = 1;

// dummy parameters as examples for creating command line parameters -> 
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
        
        // reading the number of applications of local search 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps")==0) n_apps = atoi(argv[++iarg]); 
        
        // example for creating a command line parameter param1 -> 
        // integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) {
            dummy_integer_parameter = atoi(argv[++iarg]); 
        }
        // example for creating a command line parameter param2 -> 
        // double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) {
            dummy_double_parameter = atof(argv[++iarg]);  
        }
        iarg++;
    }
}

int heuristic (const set<int>& estat) {
    return estat.size(); 
}


bool segueix_dominant(const set<int>& set, int eliminat) {  //comprovar si graf segueix sent dominant despres deliminar un
    for (int i : neighbors[eliminat]) { 


        if (((veins_dominants[i]-1)*1.0) < ceil((neighbors[i].size())/2.0) ) {
            return false; 
        }

    }

    return true; 


}


vector<int> eliminar_disponibles (const set<int>& estat)  {
    cout << "ENTRANT ELIMINA DISPONIBLES   "; 
    vector<int> retorna; 
    set<int> e = estat;
    for (int i : estat) {
        e.erase(i);
        if (segueix_dominant(e, i)) {
            retorna.push_back(i); 
        }
        e.insert(i); 
    }
    cout << "ACABAT ELIMINA DISPONIBLES   "; 
    return retorna; 
}

vector<int> add_disponibles (const set<int>& estat) {
    cout << "ENTRANT ADD DISPONIBLES   "; 
    vector<int> retorna; 
    for (int i = 0; i < neighbors.size(); ++i) {
        if (estat.find(i) == estat.end()) {
            retorna.push_back(i); 
        }
    }
    cout << "ACABANT ADD DISPONIBLES   "; 
    return retorna; 
}

set<int> genera_successor (const set<int>& estat, bool& es_eliminar, int& eliminat) {
    cout << "ENTRANT GENERA SUCCESSOR   "; 
    set<int> estat_seguent = estat; 
    vector<int> elimina = eliminar_disponibles(estat); 
    vector<int> add = add_disponibles(estat); 
    int aleatori = rand() % ((elimina.size() + add.size())-1) + 0;
    if (aleatori < elimina.size()) {
        estat_seguent.erase(elimina[aleatori]); 
        es_eliminar = true; 
        eliminat = elimina[aleatori]; 
    }
    else {
        estat_seguent.insert(add[aleatori-elimina.size()]); 
        es_eliminar = false; 
    }
    cout << "ACABAT GENERA SUCCESSOR    "; 
    return estat_seguent; 

}

double canvi_temperatura (double tempactual, double k, double lambda) {
    return (k * (pow(exp(1.0), (-lambda)*(tempactual)*1.0))  ); 
}


void decrementa_veins_dominants (int eliminat) {
    for (int i : neighbors[eliminat]) {
        veins_dominants[i]--; 
    }
}

void incrementa_veins_dominants (int afegit) {
    for (int i : neighbors[afegit]) {
        veins_dominants[i]++; 
    }
}



/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. 
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // vectors for storing the result and the computation time 
    // obtained by the <n_apps> applications of local search
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

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

    // main loop over all applications of local search
    for (int na = 0; na < n_apps; ++na) {

        // the computation time starts now
        Timer timer;

        // Example for requesting the elapsed computation time at any moment: 
        // double ct = timer.elapsed_time(Timer::VIRTUAL);

        cout << "start application " << na + 1 << endl;

        // HERE GOES YOUR LOCAL SEARCH METHOD

        // The starting solution for local search may be randomly generated, 
        // or you may incorporate your greedy heuristic in order to produce 
        // the starting solution.
        
        // Whenever you move to a new solution, first take the computation 
        // time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: 
        // cout << "value " << <value of the current solution>;
        // cout << "\ttime " << ct << endl;

        // When a local minimum is reached, store the value of the 
        // corresponding solution in vector results: 
        // results[na] = <value of the local minimum>;
        
        // Finally store the needed computation time (that is, the time 
        // measured once the local minimum is reached) in vector times: 
        // times[na] = ct;


        set<int> estat; 
        cout << "ABANS DE CREAR ESTAT INICIAL" << endl ; 
        for (int i = 0; i < neighbors.size(); i++) {
            estat.insert(i); //creacio estat inicial
            veins_dominants.push_back(neighbors[i].size()); //al principi tot els veins dominants
        }
        cout << "DESPRES DE CREAR ESTAT INICIAL" << endl; 

        int iteraciones = 4000;
        int itpertemp = 50; 
        double k = 1; 
        double lambda = 0.0000000001; 


        //AQUI TENIM ESTAT INICAL
        int contcanvitemp = 0; 
        double temp = neighbors.size();            //temperatura inicial
        for (int i = 0; i < iteraciones; i++) {
            cout << "INICI ITERACIO " << i << "     "; 
            if (contcanvitemp == itpertemp) {
                contcanvitemp = 0; 
                temp = canvi_temperatura(temp, k, lambda); 
            }
    
            bool es_eliminar = false; 
            int eliminat_afegit = 0; 
            set<int> estatseguent = genera_successor(estat, es_eliminar, eliminat_afegit); 
            double diferencia = heuristic(estat) - heuristic(estatseguent); 


            if (diferencia > 0) {
                estat = estatseguent; 
                if (es_eliminar) {
                    decrementa_veins_dominants(eliminat_afegit); 
                }
                else {
                    incrementa_veins_dominants(eliminat_afegit); 
                }
            }
            else {
                double probabilidad = pow(exp(1.0), (diferencia) / (temp * 1.0) ); 
                if ((rand() % 1000) < (int(probabilidad*1000)) ) {
                    estat = estatseguent; 
                    if (es_eliminar) {
                        decrementa_veins_dominants(eliminat_afegit); 
                    }
                    else {
                        incrementa_veins_dominants(eliminat_afegit); 
                    }
                }
            }
            contcanvitemp++; 
            cout << "FINAL ITERACIO " << i << endl; 
        }
        double ct = timer.elapsed_time(Timer::VIRTUAL); 
        cout << "VALUE: " << estat.size() <<  "\ttime " << ct << endl; 
        results[na] = estat.size(); 
        cout << endl; 
        cout << "end application " << na + 1 << endl;
    }

    // calculating the average of the results and computation times, and 
    // their standard deviations, and write them to the screen
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
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

