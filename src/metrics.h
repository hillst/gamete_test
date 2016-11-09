/*
 * PairwiseMetric:
 *     Metrics have a result object (soem data container) and a method that updates said data structure
 */ 
class PairwiseMetric{
    public:

        PairwiseMetric();
        ~PairwiseMetric();
        PairwiseMetric(const PairwiseMetric&);
        PairwiseMetricy& operator=(const PairwiseMetric&);

    public:
        virtual void compute_metric(row i, row j);

    private:
        

}
