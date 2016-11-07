#include <stdio.h>
#include <stdlib.h>
#include <boost/tokenizer.hpp>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <limits.h>

#define MIN_OBS 10
/**
 * Code for generating the data for gamete test:

    Inputs: csv file containing all sites and all states
    Outputs: counts for all pairs of state (by index? location?)

    TODO:
        - Refactor CSV reador
        - Build into a real cpp program and not this hack
        - Should be able to do the whole thing in cpp so why not
        - Rework CSVRow to be less shitty ([] override with a vector inside.... really?)

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

long int site_distance(string a, string b){
    /*
     * site_distance - accepts string a and b, computes the distance between them. 
     * 
     * assumptions - chrN-000000123456 , or chrNN-00000123456 NO scaffolds. This is RRBS only
     */
    const char * a_chr = a.c_str();
    const char * b_chr = b.c_str();
    int offset;
    if (a_chr[3] == b_chr[3] && a_chr[4] == b_chr[4]){
        if (a_chr[4] == '-'){
            offset = 5;
            const char * a_pos = a_chr + offset;
        } else{ 
            offset=6; 
        }
        long int a_pos = strtol(a_chr + offset, NULL, 10);
        if (b_chr[4] == '-'){
            offset = 5;
        } else{ 
            offset = 6; 
        }
        long int b_pos = strtol(b_chr + offset, NULL, 10);
        return abs(a_pos - b_pos);

    } else{
        return LONG_MAX;
    }
   

}
float compute_mi(CSVRow* i, CSVRow* j, int * res, int * mi){
    /*
     * Returns a count vector for the possible haplotypes (00,01,10,11)
     * PMI: H(X) − H(X|Y )
     * H(X) + Information(X), H(X|Y) is given Y (whole message)
     * so we need to see all of X and all of Y (which seems reasonable)
     * H(X, Y) is only when both X and Y are observed so there is an extra step
     * we also want to penalize missing data potentially (discount factor?)
     */
    float H_x = 0.0;
    float H_xgy = 0.0;
    int n_xy[4];
    int n_xgy[4];
    for (int k = 0; k < 4; k++){ res[k]=0; n_xy[k] =0; n_xgy[k] = 0;}
    int count = 0;
    for (int k = 1; k < i->size(); k++){ // start at 1 to skip the site-name column
        // im leaving you a puzzle
        // number of times puzzle tripped me up (to date): 1
        //
        // this loop technically works for PMI and every othe thing we've done.
        int a = ((*i)[k].c_str()[0] - 48); // 0 = 0, 1 = 1, 2 = missing
        int b = ((*j)[k].c_str()[0] - 48); // im a bitmaster
   
        
 
        if (!((b>>1) | (a>>1))){
            count++; // neither miss
            //dont we also need p given a=1 a = 0 | b=1 or b=2
            n_xgy[(a << 1) + b]++; // b = 0 or 1 so a<<1 = 2 if a == 1, 00, 01, 10, 11 is order
            n_xy[a]++; // p(x)
            n_xy[b]++; // p(y)
        } else if( !(a >> 1) ){
            // a not missing
            n_xy[a]++;
        }  else if( !(b >> 1) ){
            // b nost missing
            n_xy[b]++;
        }
    } //so now we need P(X), P(X|Y)
    float p_x0 = n_xy[0] / (n_xy[0] + n_xy[1]);
    float h_x = p_x0 * log2(p_x0) + (1-p_x0) * log2(1-p_x0);
    float p_x0_y0 = n_xgy[0] / (n_xgy[0] + n_xgy[2]); // p x = 0 given y = 0, so n(x=0,y=0)/ n(y=0,y=1)
    float p_x1_y0 = n_xgy[1] / (n_xgy[0] + n_xgy[2]); 
    float p_x0_y1 = n_xgy[0] / (n_xgy[1] + n_xgy[3]);
    float p_x1_y1 = n_xgy[1] / (n_xgy[1] + n_xgy[3]); 
    float h_xgy = p_x0_y0 * log2(p_x0_y0) + p_x1_y0 * log2(p_x1_y0) +
                  p_x0_y1 * log2(p_x0_y1) + p_x1_y1 * log2(p_x1_y1);
      
    //set return variable
    *mi = h_x - h_xgy;

    return count;
}


int compute_gamete_stats(CSVRow* i, CSVRow* j, int * res){
    /*
     * Returns a count vector for the possible haplotypes (00,01,10,11)
     * PMI: H(X) − H(X|Y )
     * H(X) + Information(X), H(X|Y) is given Y (whole message)
     * so we need to see all of X and all of Y (which seems reasonable)
     */
    for (int k = 0; k < 4; k++){ res[k]=0; }
    int count = 0;
    for (int k = 1; k < i->size(); k++){ // start at 1 to skip the site-name column
        // im leaving you a puzzle
        // number of times puzzle tripped me up (to date): 1
        //
        // this loop technically works for PMI and every othe thing we've done.
        int a = ((*i)[k].c_str()[0] - 48); // 0 = 0, 1 = 1, 2 = missing
        int b = ((*j)[k].c_str()[0] - 48); // im a bitmaster

        if (!(b>>1)){ // if b is not missing? oh right this is for p_obs_i
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
    
    auto hist_close = new int[data[0].size()][4];
    auto hist_far = new int[data[0].size()][4];
    for (int i = 0; i < data[0].size(); i++){
        for (int j = 0; j < 4; j++){
            hist_close[i][j] = 0;
            hist_far[i][j] = 0;
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
                if (site_distance(site_lookup.at(i), site_lookup.at(j)) < 100){
                    hist_close[count][n_alleles-1]++;
                } else{
                    hist_far[count][n_alleles-1]++;
                }
               // hist_close[count][n_alleles-1]++; // alleles are in 1,2,3,4 so -1
                                            // count is off by one (which is what we want)
            }
            delete [] result;
        }
     //   results.push_back(total);

        cerr << "processed outer " << i << " " << (i+1) / (time(NULL) - start + .0001) <<  " sites/second \r";
    }

    //print histogram routine 
    cout << "#hist_close" << endl;
    for (int i =0; i < data[0].size(); i++){
        for ( int j = 0; j < 4; j++){
            cout << hist_close[i][j] << " ";
        }   
        cout << endl;
    }
    cout << "#hist_far" << endl;
    for (int i =0; i < data[0].size(); i++){
        for ( int j = 0; j < 4; j++){
            cout << hist_far[i][j] << " ";
        }   
        cout << endl;
    }
    
    cerr << "\nprocessed at: " << (data.size()) / (time(NULL) - start + .0001) <<  " sites/second \n";  
    cerr << "elapsed: " << time(NULL) - start << endl;
    cerr << "elapsed: " << (time(NULL) - start) / 1000<< endl;
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
