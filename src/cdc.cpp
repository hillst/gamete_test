#include <stdio.h>

#include <string>
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <unordered_map>

#ifndef CSVROW 
    #define CSVROW 1
    #include "CSVRow.h"
    #include "metrics.h"
#endif

#define MIN_OBS 2

#ifndef GAMMA
    #define GAMMA .0000001f
#endif
/**
 * Code for generating the data for gamete test:

    Inputs: csv file containing all sites and all states
    Outputs: counts for all pairs of state (by index? location?)

    TODO:
        - Rework CSVRow to be less shitty ([] override with a vector inside.... really?)

        -OK we need a metric reader that basicaly computes stats for all the pairs, then it passes that to 

*/
using namespace std;


void load_data(char * filename, vector<CSVRow> &data, vector<string>  &site_lookup){
    ifstream file(filename); //fastest test
    cerr << "Loading data..." << endl;
    int idx = 0;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop){
        if (idx == 0){
            idx++;
        } else if ((*loop).count_present() <= MIN_OBS){
            continue;
        } else{
            data.push_back((*loop));
            site_lookup.push_back((*loop)[0]);
        }
    }
}

int main(int argc, char** argv){
        
    ifstream file(argv[1]); 
    vector<CSVRow> data;
    vector<string> site_lookup;

    // i guess we reall want distance => 


    // so we should add to distance lookup, distance => i (we have discrete distances)
    load_data(argv[1], data, site_lookup);

    unordered_map<int, int> concordant; // distancce => count concordant
    unordered_map<int, int> discordant; // distance => count discordant

    // we have i * j / 2 sites, which represent a discrete nmber of distances
    //    can we write to disk and do other stuff later? I think so.
    //

    // so what?
    //     do some stuff
    int n_sites = data.size();

    time_t start = time(NULL);
    cerr << "processing "<< data.size() << " sites" << endl;
    int n_cells = data.at(0).size(); 
    int count;

    for (int i = 0; i < data.size(); i++){
        CSVRow left = data.at(i);
        /* Some ways to make parallel perform better:
         *    no critical section
         *    iterate by num threads so they are all guaranteed to be independent (distances that is)
                - distance function and loop need to guarantee that dist(i, j) != dist(i, j+n | j+n < dist(i,j))
        */     
        for (int j = i+1; j < i+5000 && j < data.size(); j++){
            int distance = site_distance(site_lookup.at(i), site_lookup.at(j));
            if (concordant.count(distance) <= 0){
                concordant[distance] = 0;
            } 
            if (discordant.count(distance) <= 0){
                discordant[distance] = 0;            
            }
            // it's also probably bteter to do parallel work on the sites (compute distance, add to histogram)
            //    we can do them perfectly in parallel but it's worth exploring
            int cord = 0;
            int dcord = 0;
            for (int k = 1; k < n_cells; k++){
                int a = (left[k].c_str()[0] - 48);
                int b = (data[j][k].c_str()[0] - 48);
                if (((a >> 1) | (b >> 1)))  { // if missing, continue!
                    continue;
                }
                if (a == b){
                    cord++;
                } else{
                    dcord++;
                }
            }       
            // still not perfect so we should consider other approaches
            concordant[distance] += cord;
            discordant[distance] += dcord;
        }
        
        cerr << "processed outer " << i << " " << (i+1) / (time(NULL) - start + .0001) <<  " pairs/second \r";
    }
    cout << "#concordant" << endl;
    for(auto iterator = concordant.begin(); iterator != concordant.end(); iterator++) {
        cout << iterator->first << "\t" << iterator->second << endl; 
    }
    cout << "#discordant" << endl;
    for(auto iterator = discordant.begin(); iterator != discordant.end(); iterator++) {
        cout << iterator->first << "\t" << iterator->second << endl; 
    }

    //print histogram routine should go here

    return 0;
}
