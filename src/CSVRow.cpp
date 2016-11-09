#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * Some really average code for computing the row in a CSV.
 * Not super important other than the data's alignment in memory matters a lot.
 *   Can look into caching in the future.
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

