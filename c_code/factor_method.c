#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"

// number of primes below BOUND
#define NUM_PRIMES 1862
// smoothness bound
#define BOUND 16000
// only look for numbers smaller than this many bits
#define CUTOFF 257



int main(int argc, char **argv) {
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
    
    mpz_t upper_bound;
    mpz_init_set_str(upper_bound, "1", 10);
    mpz_mul_2exp(upper_bound, upper_bound, CUTOFF);

    mpz_t value; // value of p^n q that we check
    mpz_init_set_str(value, "1", 10);
    mpz_t lower; // value -1, for checking smoothness
    mpz_init(lower);
    mpz_t higher; // value +1, for checking smoothness
    mpz_init(higher);    

    for (int i = 0; i < NUM_PRIMES; i++) {
        printf("Solved %d primes\n", i);
        for (int j = 0; j < NUM_PRIMES; j++) {
            mpz_mul(value, primes[i], primes[j]);
            while (mpz_cmp(value, upper_bound) <= 0) {
                mpz_sub_ui(lower, value, 1);
                mpz_add_ui(higher, value, 1);
                if (is_smooth(lower, primes, NUM_PRIMES)) {
                    mpz_out_str(fp, 10, lower);
                    fputs("\n", fp);
                }
                if (is_smooth(higher, primes, NUM_PRIMES)) {
                    mpz_out_str(fp, 10, value);
                    fputs("\n", fp);
                }
                mpz_mul(value, value, primes[i]);
            }
        }
    }
    
    fclose(fp);

    diff_time = clock() - start_time;

    int msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    return EXIT_SUCCESS;

}

