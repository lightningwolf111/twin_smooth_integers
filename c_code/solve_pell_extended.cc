#include <gmp.h>
#include <gmpxx.h>
#include <stdbool.h>
#include <stdio.h>
#include <map>
#include "helpers.h"

using std::map;

#define SOLVED_UP_TO 2000000000

static map<mpz_class, mpz_class> solved_equations;
static bool map_initialized = false;

void solve_pell_extended(mpz_t d, mpz_t b, mpz_t result,
		mpz_t primes[], int num_primes) { 
    mpz_class d_class(d);
    if (solved_equations.count(d_class) > 0) {
        mpz_set(result, (solved_equations.at(d_class)).get_mpz_t());
	return;
    }

    solve_pell(d, b, result, primes, num_primes);
    solved_equations.insert({d_class, mpz_class(result)});
    return;

}
