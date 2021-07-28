#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

// Fills the output parameter ptr with mpz_t values of the prime numbers
// less than or equal to b. Requires that ptr point to an array with the
// space needed.
extern void primes_up_to_b(mpz_t* ptr, mpz_t b);
// Returns true if and only if arg is smooth, i.e. it factors completely
// in the set of primes given. num_primes must be the length of primes.
extern bool is_smooth(mpz_t arg, mpz_t primes[], int num_primes);
// Adds the polynomial term coeff * (eval ^ pow) into current.
extern void add_term(int coeff, int pow, mpz_t current, mpz_t eval);
// Returns true if and only if the given integer m corresponds to
// the nth solution of a Pell equation where the 4n th solution
// also gives a smooth pair.
extern bool check_fourth_poly(mpz_t m, mpz_t primes[]);
// Returns true if and only if the given integer m corresponds to
// the nth solution of a Pell equation where the 2n th solution
// also gives a smooth pair.
extern bool check_second_poly(mpz_t m, mpz_t primes[]);
// Exactly like the two above, but checks from n to 6n.
extern bool check_sixth_poly(mpz_t m, mpz_t primes[]);
// Finds smooths in the range [min, max) that factor in the array
// of primes given (who's length must be num_primes). Returns this
// through the output parameter res, which has values set to 1 for
// smooths (and unchanged values otherwise). Res should have length
// max - min, so that the ith element corresponds to the integer
// (min + i).
extern void smooths_in_range(mpz_t primes[], long min, long max,
		int num_primes, char* res);
