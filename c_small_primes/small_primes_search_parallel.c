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

// Result files
FILE *files[NUM_THREADS];

clock_t start_time, diff_time;

// The main search function, modified to put multithreading at depth n i.e.
// after n sequential recursive calls. n should be between 1 and the max number of primes
// allowed in the coefficient.
void search_n_sequential(int numFacts, mpz_t minS, mpz_t maxS, mpz_t b, mpz_t primes[], int n, int coeff_vector[], int fixed);
// Searches the space of coefficients with numFactors factors and of
// bit length between minBits and (minBits + 15).
// Puts results in the given file pointer.
// Just like the above, but leaves all computation in a single thread.
void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int start[], FILE *fp, mpz_t b, mpz_t primes[], int fixed);

int main(int argc, char **argv)
{
    int minSize, maxSize, numFacts;
    long counter_all = 0;
    clock_t solving_time_all = 0;

    printf("Minimum Coefficient size: ");
    gmp_scanf("%d", &minSize);
    printf("Maximum Coefficient size: ");
    gmp_scanf("%d", &maxSize);
    printf("Number of distinct primes in coefficient: ");
    gmp_scanf("%d", &numFacts);
    printf("Number of primes until parallelism: ");
    gmp_scanf("%ld", &num_until_par);

    start_time = clock();
    for (int i = 0; i < NUM_THREADS; i++)
    {
        solving_time[i] = 0;
        counter[i] = 0;
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
    mpz_ui_pow_ui(maxS, 2, maxSize);

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
        solving_time_all += solving_time[i];
    }
    printf("Total number of equations solved: %ld\n", counter_all);

    diff_time = clock() - start_time;
    long msec = diff_time * 1000 / CLOCKS_PER_SEC;
    printf("Total wall time: %ld s %ld ms\n", msec / 1000, msec % 1000);

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
        gmp_printf("Search called.  %d left to fix before parallelism, last prime fixed %Zd \n", n, primes[coeff_vector[fixed - 1]]);
    }
    else
    {
        printf("Search called. %d left to fix before parallelism. \n", n);
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
                coeff_vector_thread[fixed] = nextPos + 1 + thread;
                search_sequential(thread, numFacts, minS, maxS, coeff_vector_thread, files[thread], b, primes, fixed + 1);
            }
        }
    }
    else
    { // Make recursive call with n smaller.
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


void search_sequential(int thread, int numFacts, mpz_t minS, mpz_t maxS, int coeff_vector[], FILE *fp, mpz_t b, mpz_t primes[], int fixed)
{
    // int thread = omp_get_thread_num();
    // printf("[%d] Checking fixed: %d vector 0 : %d vector 1 %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
    // printf("[%d] fixed = %d, cv[0] = %d, cv[1] = %d \n", thread, fixed, coeff_vector[0], coeff_vector[1]);
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
        double start_time_solve = clock();
        // printf("[%d] solving\n", thread);
        solve_pell(current, b, result, primes, NUM_PRIMES);
        solving_time[thread] += clock() - start_time_solve;

        counter[thread]++;
        if (mpz_cmp_si(result, 0) != 0)
        {
            mpz_out_str(fp, 10, result);
            fputs(" 1\n", fp);
            // gmp_printf("Result %Zd coeff: %Zd \n", result, current);
            check_and_compute_higher_solutions(result, primes, NUM_PRIMES, fp);
            fflush(fp);
        }
        mpz_clear(result);
        if (counter[thread] % 10000 == 0)
        {
            printf("Equations solved by thread %d : %ld\n", thread, counter[thread]);

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
