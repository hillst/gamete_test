X CHROMSOME ONLY:

we want to see how well ourprogram works. Here are the steps:

1) mrbayes with n iteratoins on a naively filtered x-chromosome file, for the first pass we will filter on n_sites <= 15 (just to get something tractable and quick)

2) our program with the same parameters, generate a huge table of entries
    - goal is to interpret these tables, can we look at each site and see how they look genome wide?
    - ideally whatever policy we come up with for filtering sites from this file, we expect mrbayes to have significantly better convergence since we only include sites which consistently pass our gamete test.

    Some policies to try:
        - What is the *expected* gamete frequencies for each site?
        - Does this pass some threshold we set? (.9 maybe.. consistent in .9 pairs) (ok this is the only one i can come up with)
        - What is the distribution of gamete frequencies for each site? 
        - regions with uniform distiribution of all 4 alleles have frequent back mutation rates (altho doesnt help phyl)
        
