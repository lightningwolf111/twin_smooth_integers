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

// Exactly like the two above, but checks from n to 7n.
extern bool check_seventh_poly(mpz_t m, mpz_t primes[]);

// returns the value of the second pell solution, assuming that m
// corresponds to a fundamental solution for which the second solution
// gives a smooth pair. Result is return parameter.
extern void value_of_second_poly(mpz_t m, mpz_t result);

// returns the value of the fourth pell solution, assuming that m
// corresponds to a fundamental solution for which the fourth solution
// gives a smooth pair. Result is return parameter.
extern void value_of_fourth_poly(mpz_t m, mpz_t result);

// Finds smooths in the range [min, max) that factor in the array
// of primes given (who's length must be num_primes). Returns this
// through the output parameter res, which has values set to 1 for
// smooths (and unchanged values otherwise). Res should have length
// max - min, so that the ith element corresponds to the integer
// (min + i).
extern void smooths_in_range(mpz_t primes[], long min, long max,
		int num_primes, char* res);

// returns 0 if the pell equation with coefficient d gives no b-smooth pairs,
// and otherwise returns m where (m, m+1) are b-smooth.
// requires: d is not a square.
// result is return parameter. Must be initialized.
extern void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes);

// Returns true if and only if x * x - d * y * y == 1.
extern bool check_pell(mpz_t x, mpz_t y, mpz_t d);

// Returns through the output parameter result the squarefree part of the argument.
// Requires the argument to be smooth (within the given set of primes) and result
// to be initialized.
extern void square_free_part(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes);

// Returns through the output parameter the D from our set that gives this pair.
// Requires the argument to be smooth (within the given set of primes) and result
// to be initialized.
extern void extract_D(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes);


// The file parameter should point to a file with smooth numbers m from pairs (m, m+1)
// of smooth numbers. It creates a new file with the "_with_coeffs" appended to the 
// same name, where each line contains the coefficient of the Pell equation that gives
// this pair, followed by a space, and m.
extern void generate_file_with_coefficients(char* file_name);
