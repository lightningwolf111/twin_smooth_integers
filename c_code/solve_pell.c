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
// number of primes below 32000
#define NUM_PRIMES 3432
// smoothness bound
#define BOUND 32000

long counter;

// returns 0 if the pell equation with coefficient d gives no b-smooth pairs,
// and otherwise returns m where (m, m+1) are b-smooth.
// requires: d is not a square.
// result is return parameter. Must be initialized.
void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes);

// Returns true if and only if x * x - d * y * y == 1.
bool check_pell(mpz_t x, mpz_t y, mpz_t d);

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
        if (num_steps % 1000 == 0) {
            printf("Steps Complete: %ld\n", num_steps);
        }
    }
    
    printf("Equations solved: %ld", counter);

    fclose(fp);

    diff_time = clock() - start_time;

    int msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);    

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
                solve_pell(d, b, result, primes, NUM_PRIMES);
                //gmp_printf("%Zd\n", result);
                if (mpz_cmp_si(result,0) != 0) {
                    mpz_out_str(fp, 10, result);
                    fputs("\n", fp);
                }
            }
        }
    }
}


void two_Qprime_in_range(long min, long max, char* res, mpz_t primes[], int num_primes) { // including min, excluding max
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

}


void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes) {
    // Uses the method of continued fractions, and the recurence relations on
    // page 382 of Rosen's book Elementary Number Theory to generate the convergents.
    // Cuts off when the numerator gets too high to save time.
    
    // lists:
    mpz_t numerators[MAX_PERIOD];
    mpz_t denominators[MAX_PERIOD];
    mpz_t p_k[MAX_PERIOD];
    mpz_t q_k[MAX_PERIOD];
    mpz_t a_k[MAX_PERIOD]; // convergents to the square root of d
    
    // mpz constants:
    mpz_t one;
    mpz_init_set_str(one, "1", 10);
    mpz_t zero;
    mpz_init_set_str(zero, "0", 10);
    mpz_t cutoff;
    mpz_init_set_str(cutoff, "1", 10);
    mpz_mul_2exp(cutoff, cutoff, BIT_CUTOFF);

    // intialization:
    mpz_set(p_k[0], zero);
    mpz_set(q_k[0], one);
    mpz_sqrt(a_k[0], d);
    mpz_set(numerators[0], zero);
    mpz_set(numerators[1], one);
    mpz_set(numerators[2], a_k[0]);
    mpz_set(denominators[0], one);
    mpz_set(denominators[1], zero);
    mpz_set(denominators[2], one);
    int index = 2;

    // main loop
    while(!check_pell(numerators[index], denominators[index], d)) {
        if (mpz_cmp(numerators[index], cutoff) >= 0) {
            mpz_set(result, zero);
            return; // this pell equation does not give useful smooths
        }
        // generate the next convergent:
        // p_(k+1)
        // set_p_k_next(p_k[index-1], a_k[index-2], q_k[index-2], p_k[index-2]);
        mpz_mul(p_k[index-1], a_k[index-2], q_k[index-2]);
        mpz_sub(p_k[index-1], p_k[index-1], p_k[index-2]);
        // q_(k+1)
        // set_q_k_next(q_k[index-1], d, q_k[index-2], p_k[index-1]);
        mpz_t psquared;
        mpz_init(psquared);
        mpz_mul(psquared, p_k[index-1], p_k[index-1]);
        mpz_sub(q_k[index-1], d, psquared);
        mpz_divexact(q_k[index-1], q_k[index-1], q_k[index-2]);
        // a_(k+1)
        // set_a_k_next(d, p_k[index-1], q_k[index-1]);
        // use the fact that floor(sqrt(d)) is a_k[0].
        mpz_t p_plus_root_d;
        mpz_init(p_plus_root_d);
        mpz_add(p_plus_root_d, p_k[index-1], a_k[0]);
        mpz_fdiv_q(a_k[index-1], p_plus_root_d, q_k[index-1]);

        // generate new numerator and denominator:
        mpz_mul(numerators[index+1], numerators[index], a_k[index-1]);
        mpz_add(numerators[index+1], numerators[index+1], numerators[index-1]);

        mpz_mul(denominators[index+1], denominators[index], a_k[index-1]);
        mpz_add(denominators[index+1], denominators[index+1], denominators[index-1]);
        
        index++;
    }
    if (is_smooth(denominators[index], primes, num_primes)) {
        mpz_sub(result, numerators[index], one);
        mpz_divexact_ui(result, result, 2);
        return;
    }
    mpz_set(result, zero);
    return; // y was not smooth
}

bool check_pell(mpz_t x, mpz_t y, mpz_t d) {
    mpz_t xsquared;
    mpz_init(xsquared);
    mpz_mul(xsquared, x, x);
    mpz_t dysquared;
    mpz_init(dysquared);
    mpz_mul(dysquared, y, y);
    mpz_mul(dysquared, dysquared, d);
    mpz_sub(xsquared, xsquared, dysquared);
    if (mpz_cmp_ui(xsquared, 1) == 0) {
        return true;
    }
    return false;
}
