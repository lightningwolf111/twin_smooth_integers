#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"
#include "solve_pell_extended.h"


// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768


/*
 * Recursively checks all coefficients that can be reached my successively multiplying
 * primes into the current coefficient such that each coefficient along the way has
 * a corresponding pell equation which gives a pair of smooth numbers.
 */
void check_from(mpz_t current, mpz_t primes[], FILE *fp, mpz_t b);


int main(int argc, char **argv) {
    char* start;
    printf("Starting coefficient: \n");
    gmp_scanf("%s", start);

    clock_t start_time = clock(), diff_time;

    FILE *fp;

    fp = fopen("/tmp/res.txt", "w");
    if (fp == NULL) {
        perror("Couldn't open file.");
        return EXIT_FAILURE;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    mpz_t current_coefficient;
    mpz_init_set_str(current_coefficient, start, 10);

    check_from(current_coefficient, primes, fp, b);

    fclose(fp);


    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;
}

void check_from(mpz_t current, mpz_t primes[], FILE *fp, mpz_t b) {
    //gmp_printf("Current is: %Zd\n", current);
    mpz_t newCoeff;
    mpz_init(newCoeff);
    mpz_t result;
    mpz_init(result);

    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_mul_si(newCoeff, primes[i], 2);
        // check that coefficient is in 2QPrime:
        if (mpz_divisible_p(current, newCoeff) != 0) {
            continue;
	}
        mpz_mul(newCoeff, current, primes[i]);
	//gmp_printf("Testing new: %Zd\n", newCoeff);
	if (mpz_perfect_square_p(newCoeff) == 0) {
	    solve_pell_extended(newCoeff, b, result, primes, NUM_PRIMES);
	    //gmp_printf("%Zd\n", result);
            if (mpz_cmp_si(result,0) != 0) {
                mpz_out_str(fp, 10, result);
                fputs("\n", fp);
                check_from(newCoeff, primes, fp, b);
            }
	}
    }
    gmp_printf("Returning from current: %Zd\n", current);

    mpz_clear(newCoeff);
    mpz_clear(result);
}
