#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>
#include "config.h"

// Fills the output parameter ptr with mpz_t values of the prime numbers
// less than or equal to b. Requires that ptr point to an array with the
// space needed.
void primes_up_to_b(mpz_t* ptr, mpz_t b);

// Returns true if and only if arg is smooth, i.e. it factors completely
// in the set of primes given. num_primes must be the length of primes.
bool is_smooth(mpz_t arg, mpz_t primes[], int num_primes);

// Returns true if and only if the given integer m corresponds to
// the nth solution of a Pell equation where the 4n th solution
// also gives a smooth pair.
bool check_fourth_poly(mpz_t m, mpz_t primes[]);

// Returns true if and only if the given integer m corresponds to
// the nth solution of a Pell equation where the 2n th solution
// also gives a smooth pair.
bool check_second_poly(mpz_t m, mpz_t primes[]);

// Exactly like the two above, but checks from n to 6n.
bool check_sixth_poly(mpz_t m, mpz_t primes[]);

// Exactly like the two above, but checks from n to 7n.
bool check_seventh_poly(mpz_t m, mpz_t primes[]);

// returns the value of the second pell solution, assuming that m
// corresponds to a fundamental solution for which the second solution
// gives a smooth pair. Result is return parameter.
void value_of_second_poly(mpz_t m, mpz_t result);

// returns the value of the fourth pell solution, assuming that m
// corresponds to a fundamental solution for which the fourth solution
// gives a smooth pair. Result is return parameter.
void value_of_fourth_poly(mpz_t m, mpz_t result);

// Finds smooths in the range [min, max) that factor in the array
// of primes given (who's length must be num_primes). Returns this
// through the output parameter res, which has values set to 1 for
// smooths (and unchanged values otherwise). Res should have length
// max - min, so that the ith element corresponds to the integer
// (min + i).
void smooths_in_range(mpz_t primes[], unsigned long long min, unsigned long long max,
		int num_primes, char* res);

// Returns true if and only if x * x - d * y * y == 1.
bool check_pell(mpz_t x, mpz_t y, mpz_t d);

// Returns through the output parameter result the squarefree part of the argument.
// Requires the argument to be smooth (within the given set of primes) and result
// to be initialized.
void square_free_part(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes);

// Returns through the output parameter the D from our set that gives this pair.
// Requires the argument to be smooth (within the given set of primes) and result
// to be initialized.
void extract_D(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes);


// The file parameter should point to a file with smooth numbers m from pairs (m, m+1)
// of smooth numbers. It creates a new file with the "_with_coeffs" appended to the 
// same name, where each line contains the coefficient of the Pell equation that gives
// this pair, followed by a space, and m.
void generate_file_with_coefficients(char* file_name);

// On input m for a smooth pair (m, m+1), check all higher solutions whether they 
// correspond to a smooth pair for n = 2,3,4,...,12.
void check_and_compute_higher_solutions(mpz_t m, mpz_t primes[], int num_primes, FILE *outputfile);
void check_higher_solutions(mpz_t m, mpz_t primes[], int num_primes, bool is_pair[]);
void check_higher_solutions_upto_6(mpz_t m, mpz_t primes[], int num_primes, bool is_pair[]);

