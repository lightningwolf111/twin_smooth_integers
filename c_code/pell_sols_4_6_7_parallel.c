#include <time.h>
#include <omp.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "helpers.h"

// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768
// number of threads
#define NUM_THREADS 4

// Sieving time 
double sieving_time[NUM_THREADS];
// Time to check the 4th, 6th or 7th solution
double checking_time[NUM_THREADS];

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
	
	// optimize finding twin smooths	
	
	FILE *fp[NUM_THREADS];

	long cuts[NUM_THREADS + 1];
    long sub_interval = (end - start + (NUM_THREADS-1))/NUM_THREADS;
    cuts[0] = start;
    cuts[NUM_THREADS] = end;
    for (int i = 1; i < NUM_THREADS; i++)
    {
        cuts[i] = start + i*sub_interval;
    }

	mpz_t primes[NUM_PRIMES];
	mpz_t b;
	mpz_init_set_si(b, BOUND);
	primes_up_to_b(primes, b);

	double wall_start_time = omp_get_wtime();
	
	
	#pragma omp parallel num_threads(NUM_THREADS)
    {
        int thread = omp_get_thread_num();
        sieving_time[thread] = 0;
        checking_time[thread] = 0;

        char file_path[100];
        sprintf(file_path, "results/res_pell_sols_%d_thr_%d_%ld-%ld.txt", which, thread, start, end);

        fp[thread] = fopen(file_path, "w");
        if (fp[thread] == NULL) {
            perror("Couldn't open file.");
            //return EXIT_FAILURE;
        }

		long curr_start = cuts[thread];
		while (curr_start < cuts[thread+1]) {
			char res[step+1]; // cover last pair to overlap
			for (long i = 0; i < step + 1; i++) {
				res[i] = 0;
			}
			double sieve_start_time = omp_get_wtime();
			smooths_in_range(primes, curr_start, curr_start + step + 1, NUM_PRIMES, res);
			sieving_time[thread] += omp_get_wtime() - sieve_start_time;

			double check_start_time = omp_get_wtime();
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
									fputs(str, fp[thread]);
					}
					mpz_clear(m);
				}
			}
			checking_time[thread] += omp_get_wtime() - check_start_time;
			
			curr_start += step;
			printf("[%d] Checked up to %ld.\n", thread, curr_start);
		}
		
		printf("Sieving Time (thread %d): %lf seconds\n", thread, sieving_time[thread]);
        printf("Checking Time (thread %d): %lf seconds\n", thread, checking_time[thread]);

		fclose(fp[thread]);
	}

	printf("Total Time taken %lf seconds\n", omp_get_wtime() - wall_start_time); 

	for (int i = 0; i < NUM_PRIMES; i++) {
		mpz_clear(primes[i]);
	}
	mpz_clear(b);
	
	return EXIT_SUCCESS;
}

