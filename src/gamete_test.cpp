#include <stdio.h>
#include <boost/tokenizer.hpp>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

#define MIN_OBS 4
//oof thats like 10%!
/**
 * Code for generating the data for gamete test:

    Inputs: csv file containing all sites and all states
    Outputs: counts for all pairs of state (by index? location?)

    TODO:
        - Refactor CSV reador
        - Build into a real cpp program and not this hack
        - Should be able to do the whole thing in cpp so why not
        - Rework CSVRow to be less shitty ([] override with a vector inside.... really?)

    FastTrack TODO:
        - Use fixed size pointer array instead of vector to store data
        - Used fixed size 3d array for our data collection, do everything by index.
 */
using namespace std;
class CSVRow
{
    public:
        string const& operator[](size_t index) const{
            return m_data[index];
        }
        size_t size() const{
            return m_data.size();
        }
        void readNextRow(istream& str){
            string line;
            getline(str, line);

            stringstream lineStream(line);
            string cell;

            m_data.clear();
            while(getline(lineStream, cell, ',')){
                m_data.push_back(cell);
            }
        }
        int get_total_sites() const{
            return m_data.size();
        }
        int count_missing() const{
            int count = 0;
            for (int i = 0; i < this->size(); i++){
                if (m_data[i] == "2"){
                    count += 1;
                }
            }
            return count;
        }
        int count_present() const{
            return this->get_total_sites() - this->count_missing();
        }

    private:
        vector<string> m_data;
};

istream& operator>>(istream& str, CSVRow& data){
    data.readNextRow(str);
    return str;
}  

class CSVIterator{   
    public:
        typedef input_iterator_tag iterator_category;
        typedef CSVRow  value_type;
        typedef size_t  difference_type;
        typedef CSVRow* pointer;
        typedef CSVRow& reference;

        CSVIterator(istream& str):m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator():m_str(NULL) {}

        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}

    private:
        istream* m_str;
        CSVRow m_row;
};

int compute_gamete_stats(CSVRow* i, CSVRow* j, int * res){
    /*
     * Returns a count vector for the possible haplotypes (00,01,10,11)
     */
    for (int k = 0; k < 4; k++){ res[k]=0; }
    int count = 0;
    for (int k = 1; k < i->size(); k++){ // start at 1 to skip the site-name column
        // im leaving you a puzzle
        int a = ((*i)[k].c_str()[0] - 48); // 0 = 0, 1 = 1, 2 = missing
        int b = ((*j)[k].c_str()[0] - 48); // im a bitmaster
        if (!(b>>1)){
            count++;
        }
        if (!((a >> 1) | (b >> 1))){
            res[(a << 1) + b]++;
        }
    }
    return count;

}

void print_results(vector<float *> vals) {
    for (int i =0 ; i < vals.size(); i++){
        for(int j = 0; j < 4; j++){
            cout << vals.at(i)[j] << " ";
        }
        cout << endl;
    }
}

void print_aggregate(vector<float *> vals, vector<string> site_names){
    for (int i =0 ; i < vals.size(); i++){
        cout << site_names.at(i) << " ";
        for(int j = 0; j < 4; j++){
            cout << vals.at(i)[j] << " ";
        }
        cout << endl;
    }
}

void print_hist(int ** hist, int size_i){
    for (int i = 0; i < size_i; i++){
        for(int j =0;j<4;j++){
            cout << hist[i][j] << " ";
        }
        cout << endl;
    }

}

int main(int argc, char** argv){
    ifstream file(argv[1]); //fastest test
    vector<CSVRow> data;
    cerr << "Loading data..." << endl;
    int idx = 0;
    vector<string> site_lookup;
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
    int * result = new int[4];
    // we can have a bunch of functions for our prepare and our filtering pipeline

    vector<int*> results;    
    time_t start = time(NULL);
    // we are still looking to get faster, we need speedup from somewhere of about 1 order of magnitude
    //  maybe we can get it from unrolling the loop? are i and j continugous in memory?
    cerr << "processing "<< data.size() << " sites" << endl;
    int valid_pairs = 0;
    int count;
    
    auto hist = new int[data[0].size()][4];
    for (int i = 0; i < data[0].size(); i++){
        for (int j = 0; j < 4; j++){
            hist[i][j] = 0;
        }
    }


    for (int i = 0; i < data.size(); i++){
        CSVRow left = data.at(i);
        //float * aggregate = new float[4];

        //working on this mess -- modularize this.
        //for (int k = 0; k < 4; k++){ aggregate[k] = 0; }
        
        //float zz = 0; float zo = 0; float oz = 0; float oo =0;
        #pragma omp parallel for reduction(+:count)
        for (int j = 0; j < data.size(); j++){
            float total = 0;
            int * result = new int[4];
            int count = compute_gamete_stats(&left, &data.at(j), result); // count is n_miss
            //we could be missing cache here
            int n_alleles = 0;
            for (int k = 0; k < 4; k++){
                if (result[k] > 0){
                    n_alleles += 1;
                }
                total += result[k]; // n_observations
            }

            if (total >= 4){ // min 4 pairs is our threshold 
                hist[count][n_alleles-1]++; // alleles are in 1,2,3,4 so -1
                                            // count is off by one (which is what we want)
            }
            delete [] result;
            // now they are all the same so we should probably try just counting first
        }
     //   results.push_back(total);

        cerr << "processed outer " << i << " " << (i+1) / (time(NULL) - start + .0001) <<  " sites/second \r";
    }

    //print histogram routine 
    for (int i =0; i < data[0].size(); i++){
        for ( int j = 0; j < 4; j++){
            cout << hist[i][j] << " ";
        }   
        cout << endl;
    }
    
    cerr << "\nprocessed at: " << (data.size()) / (time(NULL) - start + .0001) <<  " sites/second \n";  
    cerr << "elapsed: " << time(NULL) - start << endl;
    cerr << "elapsed: " << (time(NULL) - start) / 1000<< endl;
    cout << "count: " << count << endl;
    //print_aggregate(results, site_lookup);
}


int computation_test(vector<CSVRow> * data){
    time_t start = time(NULL);
    for (int i = 0; i < data->size(); i++){
        for (int j = 0; j < data->size(); j++){
            int dumb = 0;
            for (int k = 0; k < 109; k++){
                dumb++;
            }
        }
        cerr << "processed outer " << i << " " << (i+1) / (time(NULL) - start + .0001) <<  " sites/second \r";
    }
    return 0;
}
