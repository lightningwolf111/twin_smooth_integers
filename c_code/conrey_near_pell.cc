#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"
#include <string.h>
#include "solve_pell_extended.h"
#include <vector>

using std::vector;
using std::pair;

// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768



int main(int argc, char **argv) {
    char filename[100];
    gmp_printf("Filename of smooth pairs: \n");
    scanf("%s", filename);

    
    char converted_file[100];
    strncpy(converted_file, filename, strlen(filename) - 4);
    strcat(converted_file, "_with_coeffs.txt");

    generate_file_with_coefficients(filename);
    
    FILE *fp;

    fp = fopen(converted_file, "r");
    if (fp == NULL) {
        perror("Couldn't open file.");
        return EXIT_FAILURE;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    vector<pair<mpz_class, mpz_class>> initial_set;
    
    
    fclose(fp);

    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;
}

