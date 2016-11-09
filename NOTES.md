Gamete test is as follows:

For every pair of sites:
    - How many combinations of gametes do we see? that is:
        00 01
        10 11
    - If you see all 4 gametes then there is no perfect phylogeny.
Idea:
    - Do this for all sites and see how this looks genome wide.


Problems:
    - Filtering criteria? well we need atleast 4 cells, but more would be better!
    - Most sites that pass the test are homogenious. This makes sense, any homogenious site with any other site by definition is a perfect phylogeny. What if there are 2 gametes? that does not seem resolvable. (00 11 00),  (???)
    - Still too slow to do genome wide
    - How do we efficiently perform this gamete test?
        - Randomization?
        - Vectorization?
        - Filtering the number of sites before getting counts?

    - Does the gamete test actually help us?
        - can produce histograms to look at # shared pairs
        - can produce histograms of near vs distal
        - can produce p_fail | n_obs_i
        - can produce n_obs_ij | n_obs_i 
        - need to look at mutual information as well
        - MI v p_fail
        - MI v n_obs_i

Thoughts:
    - Is the X-chromosome alone enough sites to build a phylogeny -- if not how do we evaluate or scale?

Progress:
    - Solved segfault issues, were related to the computation of mi
    - have not addressed reduction issues
    
    



