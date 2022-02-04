#include <time.h>
#include <stdio.h>
#include <omp.h>
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
// number of threads
#define NUM_THREADS 4

// The number of Pell equations solved.
long counter[NUM_THREADS];

// Solving time (for solving pell equations)
clock_t solving_time[NUM_THREADS];

// Number of pell solutions in desired range, without regard to smoothness of y.
long numInRange[NUM_THREADS];

// Max number to solve
long numPellToSolve;

// Result files
FILE* files[NUM_THREADS];

clock_t start_time, diff_time;

// The main search function.
// Splits the search among NUMTHREADS threads.
void search(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[]);
// Searches the space of coefficients with numFactors factors and of
// bit length between minBits and (minBits + 15).
// Puts results in the given file pointer.
// Just like the above, but leaves all computation in a single thread.
void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int start[], FILE *fp, mpz_t b, mpz_t primes[], int fixed);

int main(int argc, char **argv) {
    int minSize, numFacts;
    long counter_all = 0;
    long numInRange_all = 0;
    clock_t solving_time_all = 0;

    printf("Minimum Coefficient size: \n");
    gmp_scanf("%d", &minSize);
    printf("Number of distinct primes in coefficient: \n");
    gmp_scanf("%d", &numFacts);
    printf("Cutoff for number of Pell equations to solved: \n");
    gmp_scanf("%ld", &numPellToSolve);


    start_time = clock();
    for (int i = 0; i < NUM_THREADS; i++) {
        solving_time[i] = 0;
	    counter[i] = 0;
	    numInRange[i] = 0;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    ////////////////////////////// Core algorithm ///////////////////////////////
    
    for (int i = 0; i < NUM_THREADS; i++) {
        char file_path[100];
        sprintf(file_path, "/tmp/res_%d.txt", i);
        files[i] = fopen(file_path, "w");
        if (files[i] == NULL) {
            perror("Couldn't open file.");
            exit(EXIT_FAILURE);
    	}
    }

    mpz_t minS;
    mpz_init(minS);
    mpz_ui_pow_ui(minS, 2, minSize);
    mpz_t maxS;
    mpz_init(maxS);
    mpz_ui_pow_ui(maxS, 2, minSize + 15);
    
    search(numFacts, minS, maxS, b, primes);

    for (int i = 0; i < NUM_THREADS; i++) {
        fclose(files[i]);
    }

    ///////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < NUM_THREADS; i++) {
        counter_all += counter[i];
        numInRange_all += numInRange[i];
        solving_time_all += solving_time[i];
    }
    printf("Total number of equations solved: %ld\n", counter_all);
    printf("Total number of results in range: %ld\n", numInRange_all);

    diff_time = clock() - start_time;
    long msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Total wall time: %ld seconds %ld milliseconds\n", msec/1000, msec%1000);

    long msec_solve = solving_time_all * 1000 / CLOCKS_PER_SEC;
    printf("Total solving time: %ld s %ld ms\n", msec_solve/1000, msec_solve%1000);
    long msec_solve_per_thread = msec_solve / NUM_THREADS;
    printf("Total solving time per thread: %ld s %ld ms\n", msec_solve_per_thread/1000, msec_solve_per_thread%1000);
    
    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);
    mpz_clear(minS);
    mpz_clear(maxS);

    return EXIT_SUCCESS;
}

void search(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[]) {

    printf("Search called\n");

    int coeff_vector[numFacts];
    for (int i = 0; i < numFacts; i++) {
        coeff_vector[i] = i;
    }

    for (int firstPos = 0; firstPos < NUM_PRIMES; firstPos++) {
        coeff_vector[0] = firstPos;
        #pragma omp parallel num_threads(NUM_THREADS) 
        {
            int thread = omp_get_thread_num();
            int coeff_vector_thread[numFacts];
            for (int i = 0; i < numFacts; i++) {
                coeff_vector_thread[i] = coeff_vector[i];
            }

            for (int nextPos = 1; nextPos < NUM_PRIMES; nextPos += NUM_THREADS) { 
                coeff_vector_thread[1] = nextPos + thread;
                search_sequential(thread, numFacts, minS, maxS, coeff_vector_thread, files[thread], b, primes, 1);
            }
        }
    }
}


void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int coeff_vector[], FILE *fp, mpz_t b, mpz_t primes[], int fixed) {
    //int thread = omp_get_thread_num();
    //printf("[%d] Checking fixed: %d vector 0 : %d vector 1 %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
    //printf("[%d] fixed = %d, cv[0] = %d, cv[1] = %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
    int count = 0;
    for (int i = 0; i < NUM_THREADS; i++)
        count += counter[i];
    if (count >= numPellToSolve) {
        return;
    } 
    // Check that we are not too large
    mpz_t current;
    mpz_init_set_si(current, 2);
    for (int i = 0; i < fixed; i++) {
        if (coeff_vector[i] < NUM_PRIMES) {
            mpz_mul(current, current, primes[coeff_vector[i]]);
        }
    }
    //gmp_printf("Current %Zd \n", current);

    if (mpz_cmp(current, maxS) > 0) {
        mpz_clear(current);
        return; // Coefficient too large
    }

    // Check the current vector
    if (fixed == numFacts && mpz_cmp(current, minS) > 0) {
        mpz_t result;
        mpz_init(result);
        //double start_time_solve = clock();
        //printf("[%d] solving\n", thread);
        solve_pell(current, b, result, primes, NUM_PRIMES);
        //solving_time[thread] += clock() - start_time_solve;
        counter[thread]++;
        if (mpz_cmp_si(result,0) != 0) {
            mpz_out_str(fp, 10, result);
            fputs("\n", fp);
            //gmp_printf("Result %Zd coeff: %Zd \n", result, current);
        }
        mpz_clear(result);
        if (counter[thread] * 100 % numPellToSolve == 0) {
            printf("Thread %d finished solving: %ld\n", thread, counter[thread]);
        }
        if (counter[thread] == numPellToSolve) {
            printf("Equations solved by thread %d : %ld\n", thread, counter[thread]);
            printf("Results in range: %ld\n", numInRange[thread]);

            // diff_time = clock() - start_time;
            // int msec = diff_time * 1000 / CLOCKS_PER_SEC;
            // printf("[%d] Total Time taken %d seconds %d milliseconds\n", thread, msec/1000, msec%1000);

            // int msec_solve = solving_time[thread] * 1000 / CLOCKS_PER_SEC;
            // printf("[%d] Solving Time taken %d seconds %d milliseconds\n", thread, msec_solve/1000, msec_solve%1000);
            
            return;
        }
    }
    mpz_clear(current);

    // Recursively check other vectors
    if (fixed < numFacts) {
        int nextPos = 0;
        if (fixed != 0) {
            nextPos = coeff_vector[fixed - 1] + 1;
        }
        for (;nextPos < NUM_PRIMES; nextPos++) {
            coeff_vector[fixed] = nextPos;
            search_sequential(thread, numFacts, minS, maxS, coeff_vector, fp, b, primes, fixed + 1);
        }
    }
}
