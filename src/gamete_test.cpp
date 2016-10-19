#include <stdio.h>
#include <boost/tokenizer.hpp>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
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

int* compute_gamete_stats(CSVRow* i, CSVRow* j){
    int * res = new int[4];
    for (int k = 0; k < 4; k++){ res[k]=0; }

    for (int k = 1; k < i->size(); k++){ // start at 1 to skip the site-name column
        // im leaving you a puzzle
        int a = ((*i)[k].c_str()[0] - 48); // 0 = 0, 1 = 1, 2 = missing
        int b = ((*j)[k].c_str()[0] - 48); // im a bitmaster
        
        if (!((a >> 1) | (b >> 1))){
            res[a + b]++;
        }
        /*
        if ((*i)[k] == "-" || (*j)[k] == "-"){
            continue;
        } else if((*i)[k] == "0" && (*j)[k] == "0"){
            res[0]++;
        } else if((*i)[k] == "0" && (*j)[k] == "1"){
            res[1]++;
        } else if((*i)[k] == "1" && (*j)[k] == "0"){
            res[2]++;
        } else if ((*i)[k] == "1" && (*j)[k] == "1"){
            res[3]++;
        }
        */
    }
    return res;

}

void print_results(vector<int *> vals) {
    for (int i =0 ; i < vals.size(); i++){
        for(int j = 0; j < 4; j++){
            cout << vals.at(i)[j] << " ";
        }
        cout << endl;
    }
    
}
int main(){
//    ifstream file("shm_normal_B.csv"); // full data
    //ifstream file("shm_normal_B_100k.csv"); //slowest test

    ifstream file("plop.csv"); //fastest test
    vector<CSVRow> data;
    cout << "Loading data..." << endl;
    int idx = 0;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop){
        if (idx == 0){
            idx++;
        } else if ((*loop).count_present() < 4){
            continue;
        } else{
            data.push_back((*loop));
        }

    }

    vector<int*> results;    
    time_t start = time(NULL);
    for (int i = 0; i < data.size(); i++){
        for (int j = 0; j < i; j++){
            int genotypes[4];
            int * result = compute_gamete_stats(&data.at(i), &data.at(j));
            results.push_back(result);
        }
        cerr << "processed outer " << i << " " << (i+1) / (time(NULL) - start + .0001) <<  " sites/second \r";
    }
    print_results(results);
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
