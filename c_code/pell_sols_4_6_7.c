#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "helpers.h"

// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768


int main(int argc, char **argv) {
	int which;
	long start, end, step;
	gmp_printf("Which successive solution to check (one of {4, 6, 7}): \n");
	gmp_scanf("%d", &which);
	gmp_printf("Start: \n");
	gmp_scanf("%ld", &start);
	gmp_printf("End: \n");
	gmp_scanf("%ld", &end);
	gmp_printf("Step: \n");
    gmp_scanf("%ld", &step);
	

	mpz_t primes[NUM_PRIMES];
	mpz_t b;
	mpz_init_set_si(b, BOUND);
	primes_up_to_b(primes, b);

	// optimize finding twin smooths	
	
	FILE *fp;

	fp = fopen("/tmp/res.txt", "w");
	
	long curr_start = start;
	while (curr_start < end) {
		char res[step+1]; // cover last pair to overlap
		for (long i = 0; i < step + 1; i++) {
			res[i] = 0;
		}
		
		smooths_in_range(primes, curr_start, curr_start + step + 1, NUM_PRIMES, res);
		for (long i = 0; i < step; i++) {
			if ((res[i] == 1) && (res[i+1] == 1)) {
				// printf("smooth pair: %ld \n", i + curr_start);
				mpz_t m;
				mpz_init_set_si(m, i + curr_start);
				bool found = false;
				if (which == 4) {
				    found = check_fourth_poly(m, primes);
				} else if (which == 6) {
				    found = check_sixth_poly(m, primes);
				} else if (which == 7) {
				    found = check_seventh_poly(m, primes);
				}
				
				if (found) {
                	                //printf("SMOOTH: %ld \n ", i + curr_start);
					char str[257];
					sprintf(str, "%ld\n", i+curr_start);
					// str[257] = '\n';
        	                	fputs(str, fp);
				}
				mpz_clear(m);
			}
		}
		
		curr_start += step;
		printf("Checked up to %ld.\n", curr_start);
	}
	
	fclose(fp);
	mpz_clear(b);
	for (int i = 0; i < NUM_PRIMES; i++) {
            mpz_clear(primes[i]);
	}
	return EXIT_SUCCESS;
}

