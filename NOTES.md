Gamete test is as follows:


for every pair of sites:
    How many combinations of gametes do we see? that is:
        00
        01
        10
        11

    If you see all 4 gametes then there is no perfect phylogeny.

Idea:
    do this for all sites and see how this looks genome wide.

Problems:
    filtering criteria? well we need atleast 4 cells, but more would be better!
    First pass? lets look at normal b cells

Technical problems:
    Runtime is wayyyyyy too slow
        could try a randomization algorithm
        could try higher thresholD
        lower bound of 100 sites / second
            this would take ~ 7 hours to compute but i dont know what to do with the results.


TODO
    - ok so our bitwise thing didnt work we gotta adjust
    - then we can continue making sure everything works and is efficient.
