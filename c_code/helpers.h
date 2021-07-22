#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

extern void primes_up_to_b(mpz_t* ptr, mpz_t b);
extern bool is_smooth(mpz_t arg, mpz_t primes[], int num_primes);
//extern void add_term(int coeff, int pow, mpz_t current, mpz_t eval);
//extern bool check_fourth_poly(mpz_t m, mpz_t primes[]);
//extern bool check_second_poly(mpz_t m, mpz_t primes[]);
