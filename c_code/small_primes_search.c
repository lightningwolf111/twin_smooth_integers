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

// Number of pell solutions in desired range, without regard to smoothness of y.
long numInRange;

// Max number to solve
long numPellToSolve;

// Sieving time (searching for coefficients)
clock_t sieving_time;

// Solving time (for solving pell equations)
clock_t solving_time;

clock_t start_time, diff_time;


// Recursively searches the space of coefficients with numFactors factors and of
// bit length between minBits and (minBits + 15).
// Puts results in the given file pointer.
void search(int numFacts, mpz_t minS, mpz_t maxS, int start[], FILE *fp, mpz_t b, mpz_t primes[], int fixed);


int main(int argc, char **argv) {
    int minSize, numFacts;
    printf("Minimum Coefficient size: \n");
    gmp_scanf("%d", &minSize);
    printf("Number of distinct primes in coefficient: \n");
    gmp_scanf("%d", &numFacts);

    printf("Cutoff for number of Pell equations to solve: \n");
    gmp_scanf("%dil", &numPellToSolve);

    start_time = clock();
    sieving_time = solving_time = 0;

    FILE *fp;
    counter = 0;
    numInRange = 0;

    fp = fopen("/tmp/res.txt", "w");
    if (fp == NULL) {
        perror("Couldn't open file.");
        return EXIT_FAILURE;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    ////////////////////////////// Core algorithm ///////////////////////////////

    int coeff_vector[numFacts];
    for (int i = 0; i < numFacts; i++) {
        coeff_vector[i] = i;
    }
    mpz_t minS;
    mpz_init_set_si(minS, 1);
    mpz_ui_pow_ui(minS, 2, minSize);
    mpz_t maxS;
    mpz_init_set_si(maxS, 1);
    mpz_ui_pow_ui(maxS, 2, minSize + 15);

    
    search(numFacts, minS, maxS, coeff_vector, fp, b, primes, 0);


    ///////////////////////////// //////////////////////////////////////////////

    printf("Equations solved: %ld\n", counter);

    printf("Results in range: %ld\n", numInRange);

    fclose(fp);

    diff_time = clock() - start_time;
    int msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Total Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

    int msec_solve = solving_time * 1000 / CLOCKS_PER_SEC;
    printf("Solving Time taken %d seconds %d milliseconds\n", msec_solve/1000, msec_solve%1000);


    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;

}

void search(int numFacts, mpz_t minS, mpz_t maxS, int coeff_vector[], FILE *fp, mpz_t b, mpz_t primes[], int fixed) {
    //printf("Checking fixed: %d vector 0 : %d vector 1 %d \n", fixed, coeff_vector[0], coeff_vector[1]);
	
    // Check that we are not too large
    mpz_t current;
    mpz_init_set_si(current, 2);
    for (int i = 0; i < fixed; i++) {
        mpz_mul(current, current, primes[coeff_vector[i]]);
    }
    //gmp_printf("Current %Zd \n", current);

    if (mpz_cmp(current, maxS) > 0) {
        return; // Coefficient too large
    }

    // Check the current vector
    if (fixed == numFacts && mpz_cmp(current, minS) > 0) {
        mpz_t result;
        mpz_init(result);
        solve_pell(current, b, result, primes, NUM_PRIMES);
	counter++;
        if (mpz_cmp_si(result,0) != 0) {
            mpz_out_str(fp, 10, result);
            fputs("\n", fp);
	    //gmp_printf("Result %Zd coeff: %Zd \n", result, current);
        }
	mpz_clear(result);
	if (counter == numPellToSolve) {
            printf("Equations solved: %ld\n", counter);

            printf("Results in range: %ld\n", numInRange);

	    diff_time = clock() - start_time;
    	    int msec = diff_time * 1000 / CLOCKS_PER_SEC;
            printf("Total Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

	    int msec_solve = solving_time * 1000 / CLOCKS_PER_SEC;
    	    printf("Solving Time taken %d seconds %d milliseconds\n", msec_solve/1000, msec_solve%1000);
            fclose(fp);
            exit(0);
	}
    }

    // Recursively check other vectors
    if (fixed < numFacts) {
        int nextPos = 0;
        if (fixed != 0) {
            nextPos = coeff_vector[fixed - 1] + 1;
	}
        for (;nextPos < NUM_PRIMES; nextPos++) {
            coeff_vector[fixed] = nextPos;
            search(numFacts, minS, maxS, coeff_vector, fp, b, primes, fixed + 1);
        }
    }
    mpz_clear(current);
}
