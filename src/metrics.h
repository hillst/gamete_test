/*
 * PairwiseMetric:
 *     Metrics have a result object (soem data container) and a method that updates said data structure
 */ 
#ifndef CSVROW
    #include "CSVRow.h" // !!!!!!!!!!
#endif


#include <math.h>
#ifndef GAMMA
    #define GAMMA .000001f
#endif

#include <vector>
using namespace std;

class PairwiseMetric{
    // pure virtual
    public:
       // PairwiseMetric(){};
       // virtual ~PairwiseMetric(){};
        virtual void compute_metric(const CSVRow* , const CSVRow* )=0;
        virtual void print_result()=0; // prints the result
        const int n_alleles = 4;
};

class PFailGiveMutualInfo : virtual public PairwiseMetric{
    public:
        PFailGiveMutualInfo(int);
        ~PFailGiveMutualInfo();
        PFailGiveMutualInfo(const PFailGiveMutualInfo&);
        PFailGiveMutualInfo& operator=(const PFailGiveMutualInfo&);

    public:
        void print_result();
        void compute_metric(const CSVRow *, const CSVRow *);
    
    private:
        int resolution;
        vector< vector <int> > hist_close;
        vector< vector <int> > hist_far;
};


