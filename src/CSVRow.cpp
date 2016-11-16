#ifndef CSVROW
    #define CSVROW
    #include "CSVRow.h"
#endif

/**
 * Some really average code for computing the row in a CSV.
 * Not super important other than the data's alignment in memory matters a lot.
 *   Can look into caching in the future.
 */
        string const& CSVRow::operator[](size_t index) const{
            return m_data[index];
        }
        size_t CSVRow::size() const{
            return m_data.size();
        }
        void CSVRow::readNextRow(istream& str){
            string line;
            getline(str, line);

            stringstream lineStream(line);
            string cell;

            m_data.clear();
            while(getline(lineStream, cell, ',')){
                m_data.push_back(cell);
            }
        }
        int CSVRow::get_total_sites() const{
            return m_data.size();
        }
        int CSVRow::count_missing() const{
            int count = 0;
            for (int i = 0; i < this->size(); i++){
                if (m_data[i] == "2"){
                    count += 1;
                }
            }
            return count;
        }
        int CSVRow::count_present() const{
            return this->get_total_sites() - this->count_missing();
        }

istream& operator>>(istream& str, CSVRow& data){
    data.readNextRow(str);
    return str;
}  


        CSVIterator::CSVIterator(istream& str):m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator::CSVIterator():m_str(NULL) {}

        CSVIterator& CSVIterator::operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        CSVIterator CSVIterator::operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& CSVIterator::operator*()   const       {return m_row;}
        CSVRow const* CSVIterator::operator->()  const       {return &m_row;}


        bool CSVIterator::operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool CSVIterator::operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}


