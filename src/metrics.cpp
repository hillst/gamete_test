/*
 * Metrics are here to allow for a polymorphic interface for computing metrics.
 * We need to address this compute metric problem
 *    we will be doing O(n) for all metrics which is not ideal!
 *    if our metric is stateful we could get clever there, like we havei nternal data, and we are updating.

 * WHERE WE AT
 *   We did something with our 2d vectors... gotta get that working dynamically and i'd 
 *     prefer not to have a new class. maybe structs. maybe a new class. dunno. 
 *     try and figure out how much data our current guy can hold and how to adjust that.
 */ 
#include "metrics.h"
#define N_ALLELES 4

PFailGiveMutualInfo::PFailGiveMutualInfo(int const resolution = 10){
    /**
     * OBJECTIVES: We want to be able to compute and store a metric that acts on pairs of rows. 
     * The interal data structure will update itself as we call compute_metric or whatever
     *   So compute_metric is not immutable and not sensitive to making the same calls twice.
     *   The solution could be to provide a custom iterator to use with these metric classes 
     */ 
    vector< vector < int > > hist_far;
    vector< vector < int > > hist_close;
    // resolution + 1 to make sure we can handle 0 - resolution inclusive
    hist_far.resize(resolution+1, vector<int>(this->n_alleles, 0));
    hist_close.resize(resolution+1, vector<int>(this->n_alleles, 0));

    this->hist_far = hist_far;
    this->hist_close = hist_close;
    this->resolution = resolution;
}  

 
PFailGiveMutualInfo::~PFailGiveMutualInfo(){
    /// I dont know if we need to implement htis.
}


PFailGiveMutualInfo::PFailGiveMutualInfo(const PFailGiveMutualInfo&){
    // copy constructor
}


PFailGiveMutualInfo& PFailGiveMutualInfo::operator=(const PFailGiveMutualInfo&){
    // assignment operator (similar to copy, but actually 
    //     should be copy + delete old
    // a = b
    // delete a
    // copy b to a

}

void PFailGiveMutualInfo::print_result(){
    //prints data! could be more generic if we let it take the shape of the histgram
}
void PFailGiveMutualInfo::compute_metric(const CSVRow* i, const CSVRow* j){
    /*
     * Returns a count vector for the possible haplotypes (00,01,10,11)
     * PMI: H(X) − H(X|Y )
     * H(X) + Information(X), H(X|Y) is given Y (whole message)
     * so we need to see all of X and all of Y (which seems reasonable)
     * H(X, Y) is only when both X and Y are observed so there is an extra step
     * we also want to penalize missing data potentially (discount factor?)
     */

    // mutual information is a float.
    float H_x = 0.f;
    float H_xgy = 0.f;
    float n_xy[2];
    float n_xgy[4]; //these correspond to two things: the gametes we observe and the order, so 0 | 0, 0 | 1... etc
    float p_x0=0.f; float p_x1=0.f; float p_y0=0.f; float p_y1=0.f;
    float p_x0_y0=0.f; float p_x0_y1=0.f; float p_x1_y0=0.f; float p_x1_y1=0.f;
    float portion_y0 = 0.f; float portion_y1 = 0.f;
    int n_alleles=0;
    float h_x = 0.f; float h_xgy = 0.f;
    int n_obs_i = 0;


    for (int k = 1; k < i->size(); k++){
        int a = ((*i)[k].c_str()[0] - 48);
        int b =((*j)[k].c_str()[0] - 48);

        if (!(a >> 1)){ // n_obs_i
            n_obs_i++;
        }   


        // compute n(x=0), n(x=1), n(x | y) 
        if (!((a >> 1) | (b >> 1)))  {
            n_xgy[(a << 1) + b]++; //00, 01, 10, 11 is the order, also the number of alleles which is neat
        } else if(!(a>>1)){
            n_xy[a]++;
        } else if (!(b>>1)){
            n_xy[b]++;
        }
    }
    for (int k = 0; k < 4; k++){ if (n_xgy[k] > 0) n_alleles++; } 
    // compute p(x)
    if ((n_xgy[0] + n_xgy[1] == 0)){
        p_x0 = 0.f;
    } else{
        p_x0 = float(n_xgy[0] + n_xgy[1]) / (n_xgy[0] + n_xgy[1] + n_xgy[2] + n_xgy[3]);
    }
    p_x1 = 1 - p_x0;
    // compute p(y)
    if ((n_xgy[0] + n_xgy[1] == 0)){
        p_y0 = 0.f;
    } else{
        p_y0 = float(n_xgy[0] + n_xgy[2]) / (n_xgy[0] + n_xgy[1] + n_xgy[2] + n_xgy[3]);
    }
    p_y1 = 1 - p_y0;

    //compute p(x|y=0) 
    if ( n_xgy[0] == 0 ){
        p_x0_y0 = 0.f;
    } else{
        p_x0_y0 = n_xgy[0]/ (n_xgy[0] + n_xgy[2]);
    }
    p_x0_y1 = 1 - p_x0_y1;

    // compute h_x
    if (p_x0 < GAMMA || p_x1 < GAMMA){
        h_x = 0.f; 
    } else{
        h_x = p_x0 * log2(p_x0) + (1-p_x0) * log2(1 - p_x0);
        h_x = -1 * h_x;
    }
    // compute h_x | y = 0

    if (p_x0_y0 < GAMMA || p_x1_y0 < GAMMA){
        portion_y0 = 0.f;
    } else{
        portion_y0 = p_x0_y0 * log2(p_x0_y0) + p_x1_y0 * log2(p_x1_y0);
        portion_y0 *= -1 * (p_y0);
    }

    //compute h_x | y = 1
    if (p_x0_y1 < GAMMA || p_x1_y1 < GAMMA){
        portion_y1 = 0.f;
    } else{
        portion_y1 = p_x0_y1 * log2(p_x0_y1) + p_x1_y1 * log2(p_x1_y1);
        portion_y1 *= -1 * (p_y1);
    }

    h_xgy = portion_y0 + portion_y1;
    float mi = h_x - h_xgy; // oh right we have to do our clever rounding still.
    this->hist_far[int(mi * this->resolution)][ n_alleles ]++; 
    
}

int definitely_not_main(){
    printf("hello world\n");
    PFailGiveMutualInfo * m = new PFailGiveMutualInfo(100);
    m->compute_metric(0,0); // lol why wouldnt this segfault
}
