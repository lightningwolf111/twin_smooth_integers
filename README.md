The goal of this repository is to find twin smooth integers of around 256 bits with prime sum, which can be used as public parameters in variants of the B_SIDH protocol.

By definition, for a positive integer b, an integer m is b smooth iff all of the prime factors of m are at most b. Integers (m, m+1) are called twin b-smooth integers if all of the prime factors of both m and m+1 are at most b. We are looking for large pairs of integers (m, m+1), for which 2m+1 is prime, and (m, m+1) are b-smooth for as small a b as possible.

The repository is organized as follows:
- theory contains general information about the problem of finding twin smooth integers, and explains the chosen methods.
- sage_code contains python implementations (using sage) of various methods for finding twin smooth integers.
- c_code contains more efficient c implementations (using the gmp library) of the more promising methods
- pell_solution_results contains values of m for which (m, m+1) is 32000-smooth, where these values were found by running the c code to solve Pell equations.
