#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"

// The maximum period of continued fractions that we allow. 1000 should be very safe.
#define MAX_PERIOD 1000
// The cutoff for continued fractions. Will only find pairs up to this cutoff-1 bits from
// fundamental solutions.
#define BIT_CUTOFF (258)
// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768

// The number of Pell equations solved.
long counter;

// Sieving time (searching for coefficients)
clock_t sieving_time;

// Solving time (for solving pell equations)
clock_t solving_time;

// Solves pell equations in 2 * Qprime and puts m into the file for any
// smooth pairs (m, m+1) resulting from fundamental solutions to those pell
// equations.
void sieve_interval_pell(long start, long end, FILE *fp, mpz_t b, mpz_t primes[]);

// Returns a char array through the output parameter res, where are return value of 1 in
// res means that min+index is in two*QPrime (see Lehmer's definition).
void two_Qprime_in_range(long min, long max, char* res, mpz_t primes[], int num_primes);

int main(int argc, char **argv) {
    long start, end, step;
    printf("Starting coefficient: \n");
    gmp_scanf("%ld", &start);
    printf("Ending coefficient: \n");
    gmp_scanf("%ld", &end);
    printf("Step: \n");
    gmp_scanf("%ld", &step);
    
    clock_t start_time = clock(), diff_time;
    sieving_time = solving_time = 0;

    FILE *fp;
    counter = 0;

    fp = fopen("/tmp/res.txt", "w");
    if (fp == NULL) {
        perror("Couldn't open file.");
        return EXIT_FAILURE;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);
    
    long num_steps = 0;
    long curr_start = start;
    while(curr_start < end) {
        //printf(curr_start);
        sieve_interval_pell(curr_start, curr_start + step, fp, b, primes);
        curr_start += step;
        num_steps++;
        if (num_steps % 100 == 0) {
            printf("Steps Complete: %ld\n", num_steps);
        }
    }
    
    printf("Equations solved: %ld\n", counter);

    fclose(fp);

    diff_time = clock() - start_time;

    int msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Total Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);    

    int msec_sieve = sieving_time * 1000 / CLOCKS_PER_SEC;
    printf("Sieving Time taken %d seconds %d milliseconds\n", msec_sieve/1000, msec_sieve%1000);

    int msec_solve = solving_time * 1000 / CLOCKS_PER_SEC;
    printf("Solving Time taken %d seconds %d milliseconds\n", msec_solve/1000, msec_solve%1000);
    
    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;

}

void sieve_interval_pell(long start, long end, FILE *fp, mpz_t b, mpz_t primes[]) {
    char res[end - start + 1]; // cover last pair to overlap
    for (long i = 0; i < end - start + 1; i++) {
            res[i] = 0;
    }
    two_Qprime_in_range(start, end, res, primes, NUM_PRIMES);
    for (long i = start; i < end; i++) {
        if (res[i-start] == 1) {
            counter++;
            // printf("In 2QPRIME: %ld\n", (i));
            mpz_t result;
            mpz_init(result);
            mpz_t d;
            mpz_init(d);
            mpz_set_ui(d, i);
            if (mpz_perfect_square_p(d) == 0) { // not perfect square TODO Maybe replace with removing "4" from 2QPrime
		//gmp_printf("Current is: %Zd\n", d);
                //gmp_printf("B is: %Zd\n", b);
                solve_pell(d, b, result, primes, NUM_PRIMES);

                //gmp_printf("%Zd\n", result);
                if (mpz_cmp_si(result,0) != 0) {
                    mpz_out_str(fp, 10, result);
                    fputs("\n", fp);
                }
            }
            mpz_clear(result);
            mpz_clear(d);
        }
    }
}


void two_Qprime_in_range(long min, long max, char* res, mpz_t primes[], int num_primes) { // including min, excluding max
    clock_t start_time = clock();
    long nums[max - min];
    for (long i = 0; i < max - min; i++) {
        nums[i] = min + i;
    }
    for (int i = 0; i < num_primes; i++) {
        long div_by = mpz_get_ui(primes[i]);
        long mult = min/div_by + (min % div_by != 0);
        while (mult * div_by < max) {
            nums[mult * div_by - min] /= mpz_get_ui(primes[i]);
            mult += 1;
        }
    }
    for (long i = 0; i < max - min; i++) {
        //printf("start: %ld index: %ld value: %ld", min, i, nums[i]);
        if ((nums[i] == 2 || nums[i] == 1) && (i+min)%2 == 0) {
            res[i] = 1;
        }
    }
    sieving_time += clock() - start_time;
}

