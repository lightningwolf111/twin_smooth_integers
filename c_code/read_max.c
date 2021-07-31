#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include "helpers.h"

// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768

int main(int argc, char **argv) {
	char filename[100];
    gmp_printf("Filename: \n");
	scanf("%s", filename);

	mpz_t primes[NUM_PRIMES];
	mpz_t b;
	mpz_init_set_si(b, BOUND);
	primes_up_to_b(primes, b);

	FILE *read_from;

    read_from = fopen(filename, "r");
	
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    mpz_t m;
    mpz_init(m);
    mpz_t max;
    mpz_init(max);
	while ((read = getline(&line, &len, read_from)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        mpz_set_str(m, line, 10);
        if (mpz_cmp(m, max) > 0) {
            mpz_set(max, m);
        }
    }
    gmp_printf("Maximum pair: %Zd\n", max);
    
    mpz_clear(m);
    mpz_clear(max);	

	fclose(read_from);
	mpz_clear(b);
	for (int i = 0; i < NUM_PRIMES; i++) {
            mpz_clear(primes[i]);
	}
	return EXIT_SUCCESS;
}

