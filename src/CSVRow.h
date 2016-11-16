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

class CSVRow{
    public:
        string const& operator[](size_t index) const;
        size_t size() const;
        void readNextRow(istream& str);
        int get_total_sites() const;
        int count_missing() const;
        int count_present() const;

    private:
        vector<string> m_data;
};

istream& operator>>(istream& str, CSVRow& data);

class CSVIterator{
    public:
        typedef input_iterator_tag iterator_category;
        typedef CSVRow  value_type;
        typedef size_t  difference_type;
        typedef CSVRow* pointer;
        typedef CSVRow& reference;

        CSVIterator(istream& str);
        CSVIterator(); 

        CSVIterator& operator++();
        CSVIterator operator++(int);
        CSVRow const& operator*()   const;
        CSVRow const* operator->()  const;

        bool operator==(CSVIterator const& rhs); 
        bool operator!=(CSVIterator const& rhs);

    private:
        istream* m_str;
        CSVRow m_row;
};

