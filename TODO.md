1) Create a list of metrics
2) Create a generic metric class with a single virtual metric function
    - accepts four counts
    - how many values allowed per value?
    - how do we signal cases with more than one value?

3) ok bettr idea
    - we are gonna let this go on the fly while we evaluate stuff but maybe all at once
    - compute histogram on the fly (n^2) runtime m space (m is small).
        - histogram of 1 2 3 4 gametes for each pair
        - we are just producing a count of n_cells x 4 table
            n_obs x n-gametes essentially
            1,2,3,4,5,6,7 obs => 1, 2, 3 ,4 gametes
    - ultimately we are going to need to reuse these functions to compute probabilities before filtering -- it might make more sense to make them distinct functions.
