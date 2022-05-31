#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"
#include "config.h"

// The point at which to branch into parallelism.
const int num_until_par;

// The number of Pell equations solved.
long counter[NUM_THREADS];

// Solving time (for solving pell equations)
clock_t solving_time[NUM_THREADS];

// Number of pell solutions in desired range, without regard to smoothness of y.
long numInRange[NUM_THREADS];

// Max number to solve
long numPellToSolve;

// Result files
FILE *files[NUM_THREADS];

clock_t start_time, diff_time;

// The main search function, modified to put multithreading at depth n i.e.
// after n sequential recursive calls. n should be between 1 and the max number of primes
// allowed in the coefficient.
void search_n_sequential(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[], int n, int coeff_vector[], int fixed);
// The main search function.
// Splits the search among NUMTHREADS threads.
void search(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[]);
// Searches the space of coefficients with numFactors factors and of
// bit length between minBits and (minBits + 15).
// Puts results in the given file pointer.
// Just like the above, but leaves all computation in a single thread.
void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int start[], FILE *fp, mpz_t b, mpz_t primes[], int fixed);
// returns 0 if the pell equation with coefficient d gives no b-smooth pairs,
// and otherwise returns m where (m, m+1) are b-smooth.
// requires: d is not a square.
// result is return parameter. Must be initialized.
void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes, long *numInRange);


int main(int argc, char **argv)
{
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
    printf("Number of primes until parallelism: \n");
    gmp_scanf("%ld", &num_until_par);

    start_time = clock();
    for (int i = 0; i < NUM_THREADS; i++)
    {
        solving_time[i] = 0;
        counter[i] = 0;
        numInRange[i] = 0;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    ////////////////////////////// Core algorithm ///////////////////////////////

    for (int i = 0; i < NUM_THREADS; i++)
    {
        char file_path[100];
        sprintf(file_path, "results/res_small_primes_search_%d.txt", i);
        files[i] = fopen(file_path, "w");
        if (files[i] == NULL)
        {
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

    int coeff_vector[numFacts];
    for (int i = 0; i < numFacts; i++)
    {
        coeff_vector[i] = i;
    }

    search_n_sequential(numFacts, minS, maxS, b, primes, num_until_par, coeff_vector, 0);

    for (int i = 0; i < NUM_THREADS; i++)
    {
        fclose(files[i]);
    }

    ///////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < NUM_THREADS; i++)
    {
        counter_all += counter[i];
        numInRange_all += numInRange[i];
        solving_time_all += solving_time[i];
    }
    printf("Total number of equations solved: %ld\n", counter_all);
    printf("Total number of results in range: %ld\n", numInRange_all);

    diff_time = clock() - start_time;
    long msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Total wall time: %ld seconds %ld milliseconds\n", msec / 1000, msec % 1000);

    long msec_solve = solving_time_all * 1000 / CLOCKS_PER_SEC;
    printf("Total solving time: %ld s %ld ms\n", msec_solve / 1000, msec_solve % 1000);
    long msec_solve_per_thread = msec_solve / NUM_THREADS;
    printf("Total solving time per thread: %ld s %ld ms\n", msec_solve_per_thread / 1000, msec_solve_per_thread % 1000);

    for (int i = 0; i < NUM_PRIMES; i++)
    {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);
    mpz_clear(minS);
    mpz_clear(maxS);

    return EXIT_SUCCESS;
}

void search_n_sequential(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[], int n, int coeff_vector[], int fixed)
{
    if (fixed != 0)
    {
        gmp_printf("Search called. Left to fix %d, last prime fixed %Zd \n", n, primes[coeff_vector[fixed - 1]]);
    }
    else
    {
        printf("Search called. Left to fix %d. \n", n);
    }
    if (n == 0)
    { // Use parallelism
#pragma omp parallel num_threads(NUM_THREADS)
        {
            int thread = omp_get_thread_num();
            int coeff_vector_thread[numFacts];
            for (int i = 0; i < numFacts; i++)
            {
                coeff_vector_thread[i] = coeff_vector[i];
            }
            for (int nextPos = coeff_vector[fixed - 1]; nextPos < NUM_PRIMES; nextPos += NUM_THREADS)
            {
                coeff_vector_thread[fixed] = nextPos + thread;
                search_sequential(thread, numFacts, minS, maxS, coeff_vector_thread, files[thread], b, primes, fixed + 1);
            }
        }
    }
    else
    { // Make recursive call with n smaller.
        long count = 0;
        for (int i = 0; i < NUM_THREADS; i++)
            count += counter[i];
        if (count >= numPellToSolve)
        {
            return;
        }
        // Check that we are not too large
        mpz_t current;
        mpz_init_set_si(current, 2);
        for (int i = 0; i < fixed; i++)
        {
            if (coeff_vector[i] < NUM_PRIMES)
            {
                mpz_mul(current, current, primes[coeff_vector[i]]);
            }
        }
        // gmp_printf("Current %Zd \n", current);

        if (mpz_cmp(current, maxS) > 0)
        {
            mpz_clear(current);
            return; // Coefficient too large
        }
        mpz_clear(current);
        // Recursively check other vectors
        int nextPos = 0;
        if (fixed != 0)
        {
            nextPos = coeff_vector[fixed - 1] + 1;
        }
        for (; nextPos < NUM_PRIMES; nextPos++)
        {
            coeff_vector[fixed] = nextPos;
            search_n_sequential(numFacts, minS, maxS, b, primes, n - 1, coeff_vector, fixed + 1);
        }
    }
}

void search(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[])
{

    printf("Search called\n");

    int coeff_vector[numFacts];
    for (int i = 0; i < numFacts; i++)
    {
        coeff_vector[i] = i;
    }

    for (int firstPos = 0; firstPos < NUM_PRIMES; firstPos++)
    {
        coeff_vector[0] = firstPos;
#pragma omp parallel num_threads(NUM_THREADS)
        {
            int thread = omp_get_thread_num();
            int coeff_vector_thread[numFacts];
            for (int i = 0; i < numFacts; i++)
            {
                coeff_vector_thread[i] = coeff_vector[i];
            }

            for (int nextPos = 1; nextPos < NUM_PRIMES; nextPos += NUM_THREADS)
            {
                coeff_vector_thread[1] = nextPos + thread;
                search_sequential(thread, numFacts, minS, maxS, coeff_vector_thread, files[thread], b, primes, 1);
            }
        }
    }
}

void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int coeff_vector[], FILE *fp, mpz_t b, mpz_t primes[], int fixed)
{
    // int thread = omp_get_thread_num();
    // printf("[%d] Checking fixed: %d vector 0 : %d vector 1 %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
    // printf("[%d] fixed = %d, cv[0] = %d, cv[1] = %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
    long count = 0;
    for (int i = 0; i < NUM_THREADS; i++)
        count += counter[i];
    if (count >= numPellToSolve)
    {
        return;
    }
    // Check that we are not too large
    mpz_t current;
    mpz_init_set_si(current, 2);
    for (int i = 0; i < fixed; i++)
    {
        if (coeff_vector[i] < NUM_PRIMES)
        {
            mpz_mul(current, current, primes[coeff_vector[i]]);
        }
    }
    // gmp_printf("Current %Zd \n", current);

    if (mpz_cmp(current, maxS) > 0)
    {
        mpz_clear(current);
        return; // Coefficient too large
    }

    // Check the current vector
    if (fixed == numFacts && mpz_cmp(current, minS) > 0)
    {
        mpz_t result;
        mpz_init(result);
        // double start_time_solve = clock();
        // printf("[%d] solving\n", thread);
        solve_pell(current, b, result, primes, NUM_PRIMES, &numInRange[thread]);
        // solving_time[thread] += clock() - start_time_solve;

        counter[thread]++;
        if (mpz_cmp_si(result, 0) != 0)
        {
            mpz_out_str(fp, 10, result);
            fputs("\n", fp);
            // gmp_printf("Result %Zd coeff: %Zd \n", result, current);
            check_and_compute_higher_solutions(result, primes, NUM_PRIMES, fp);
        }
        mpz_clear(result);
        if (counter[thread] * 100 % numPellToSolve == 0)
        {
            printf("Thread %d finished solving: %ld\n", thread, counter[thread]);
        }
        if (counter[thread] == numPellToSolve)
        {
            printf("Equations solved by thread %d : %ld\n", thread, counter[thread]);
            printf("Results in range: %ld\n", numInRange[thread]);

            // diff_time = clock() - start_time;
            // int msec = diff_time * 1000 / CLOCKS_PER_SEC;
            // printf("[%d] Total Time taken %d seconds %d milliseconds\n", thread, msec/1000, msec%1000);

            // int msec_solve = solving_time[thread] * 1000 / CLOCKS_PER_SEC;
            // printf("[%d] Solving Time taken %d seconds %d milliseconds\n", thread, msec_solve/1000, msec_solve%1000);
            mpz_clear(current);
            return;
        }
    }
    mpz_clear(current);

    // Recursively check other vectors
    if (fixed < numFacts)
    {
        int nextPos = 0;
        if (fixed != 0)
        {
            nextPos = coeff_vector[fixed - 1] + 1;
        }
        for (; nextPos < NUM_PRIMES; nextPos++)
        {
            coeff_vector[fixed] = nextPos;
            search_sequential(thread, numFacts, minS, maxS, coeff_vector, fp, b, primes, fixed + 1);
        }
    }
}

void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes, long *numInRange) {
    // Uses the method of continued fractions, and the recurrence relations on
    // page 382 of Rosen's book Elementary Number Theory to generate the convergents.
    // Cuts off when the numerator gets too high to save time.

    clock_t start_time = clock();
    clock_t solving_time = 0;
    
    mpz_t one;
    mpz_t zero;
    mpz_t cutoff;
    mpz_init_set_str(one, "1", 10);
    mpz_init_set_str(zero, "0", 10);
    mpz_init_set_str(cutoff, "1", 10);
    mpz_mul_2exp(cutoff, cutoff, BIT_CUTOFF);

    // lists:
    mpz_t numerators[MAX_PERIOD];
    mpz_t denominators[MAX_PERIOD];
    mpz_t p_k[MAX_PERIOD];
    mpz_t q_k[MAX_PERIOD];
    mpz_t a_k[MAX_PERIOD]; // convergents to the square root of d

    // intialization:
    mpz_init_set(p_k[0], zero);
    mpz_init_set(q_k[0], one);
    mpz_init(a_k[0]);
    mpz_sqrt(a_k[0], d);
    mpz_init_set(numerators[0], zero);
    mpz_init_set(numerators[1], one);
    mpz_init_set(numerators[2], a_k[0]);
    mpz_init_set(denominators[0], one);
    mpz_init_set(denominators[1], zero);
    mpz_init_set(denominators[2], one);
    int index = 2;

    mpz_t psquared;
    mpz_init(psquared);
    mpz_t p_plus_root_d;
    mpz_init(p_plus_root_d);

    // main loop
    while(!check_pell(numerators[index], denominators[index], d)) {
        if (mpz_cmp(numerators[index], cutoff) >= 0) {
            mpz_set(result, zero);
            for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
                mpz_clear(p_k[i]);
                mpz_clear(q_k[i]);
                mpz_clear(a_k[i]);
            }
            for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
                mpz_clear(numerators[i]);
                mpz_clear(denominators[i]);
            }
	    solving_time += clock() - start_time;
            mpz_clear(psquared);
            mpz_clear(p_plus_root_d);
	    mpz_clear(one);
            mpz_clear(zero);
            mpz_clear(cutoff);
            return; // this pell equation does not give useful smooths
        }
        mpz_init(p_k[index-1]);
        mpz_init(q_k[index-1]);
        mpz_init(a_k[index-1]);
        mpz_init(numerators[index+1]);
        mpz_init(denominators[index+1]);
        // generate the next convergent:
        // p_(k+1)
        // set_p_k_next(p_k[index-1], a_k[index-2], q_k[index-2], p_k[index-2]);
        mpz_mul(p_k[index-1], a_k[index-2], q_k[index-2]);
        mpz_sub(p_k[index-1], p_k[index-1], p_k[index-2]);
        // q_(k+1)
        // set_q_k_next(q_k[index-1], d, q_k[index-2], p_k[index-1]);
        mpz_mul(psquared, p_k[index-1], p_k[index-1]);
        mpz_sub(q_k[index-1], d, psquared);
        mpz_divexact(q_k[index-1], q_k[index-1], q_k[index-2]);
        // a_(k+1)
        // set_a_k_next(d, p_k[index-1], q_k[index-1]);
        // use the fact that floor(sqrt(d)) is a_k[0].
        mpz_add(p_plus_root_d, p_k[index-1], a_k[0]);
        mpz_fdiv_q(a_k[index-1], p_plus_root_d, q_k[index-1]);

        // generate new numerator and denominator:
        mpz_mul(numerators[index+1], numerators[index], a_k[index-1]);
        mpz_add(numerators[index+1], numerators[index+1], numerators[index-1]);

        mpz_mul(denominators[index+1], denominators[index], a_k[index-1]);
        mpz_add(denominators[index+1], denominators[index+1], denominators[index-1]);

        index++;
    }

    ////////////////////
    mpz_t minBound;
    mpz_init_set_str(minBound, "1", 10);
    mpz_mul_2exp(minBound, minBound, BIT_CUTOFF - 18); // Set 2^240 bits as the min
    if (mpz_cmp(numerators[index], minBound) > 0) {
        *numInRange++;
    }
    mpz_clear(minBound);
    ////////////////////

    if (is_smooth(denominators[index], primes, num_primes)) {
        mpz_sub(result, numerators[index], one);
        mpz_divexact_ui(result, result, 2);
        solving_time += clock() - start_time;

        for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
            mpz_clear(p_k[i]);
            mpz_clear(q_k[i]);
            mpz_clear(a_k[i]);
        }
        for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
            mpz_clear(numerators[i]);
            mpz_clear(denominators[i]);
        }
        mpz_clear(psquared);
        mpz_clear(p_plus_root_d);
        mpz_clear(one);
        mpz_clear(zero);
        mpz_clear(cutoff);
        return;
    }
    mpz_set(result, zero);
    solving_time += clock() - start_time;

    for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
            mpz_clear(p_k[i]);
            mpz_clear(q_k[i]);
            mpz_clear(a_k[i]);
    }
    for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
        mpz_clear(numerators[i]);
        mpz_clear(denominators[i]);
    }
    mpz_clear(psquared);
    mpz_clear(p_plus_root_d);
    mpz_clear(one);
    mpz_clear(zero);
    mpz_clear(cutoff);
    return; // y was not smooth
}
