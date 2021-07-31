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
    FILE *write_to_2;
    FILE *write_to_4;

	write_to_2 = fopen("/tmp/res_2nd_sols.txt", "w");
    write_to_4 = fopen("/tmp/res_4th_sols.txt", "w");
    read_from = fopen(filename, "r");
	
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    mpz_t m;
    mpz_init(m);
    mpz_t second;
    mpz_t fourth;
	while ((read = getline(&line, &len, read_from)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        mpz_set_str(m, line, 10);
        //gmp_printf("m: %Zd\n", m);
        if (check_second_poly(m, primes)) {
            mpz_init(second);
            value_of_second_poly(m, second);
            mpz_out_str(write_to_2, 10, second);
            fputs("\n", write_to_2);
            mpz_clear(second);
        }
        if (check_fourth_poly(m, primes)) {
            mpz_init(fourth);
            value_of_fourth_poly(m, fourth);
            mpz_out_str(write_to_4, 10, fourth);
            fputs("\n", write_to_4);
            mpz_clear(fourth);
        }
    }
    
    mpz_clear(m);	

	fclose(read_from);
    fclose(write_to_2);
    fclose(write_to_4);
	mpz_clear(b);
	for (int i = 0; i < NUM_PRIMES; i++) {
            mpz_clear(primes[i]);
	}
	return EXIT_SUCCESS;
}

