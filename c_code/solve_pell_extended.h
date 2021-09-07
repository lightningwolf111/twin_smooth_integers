#include <gmp.h>
#include <stdbool.h>

/*
 * A wrapper over the solve_pell function defined in helpers.h which uses the fact that we
 * have already covered all coefficients under COVERED_UP_TO to speed up the computation
 * in those cases (skip solving the pell equation, and lookup the previous result).
 *
 * This is useful for long computations. 
 *
 */
extern void solve_pell_extended(mpz_t d, mpz_t b,
		mpz_t result, mpz_t primes[], int num_primes);
