#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include "helpers.h"

// Solves pell equations in 2 * Qprime and puts m into the file for any
// smooth pairs (m, m+1) resulting from fundamental solutions to those pell
// equations.
void sieve_interval_pell(long start, long end, FILE *fp, mpz_t b, mpz_t primes[], int thread, clock_t* sieving_time, clock_t* solving_time, long* counter);

// Returns a char array through the output parameter res, where are return value of 1 in
// res means that min+index is in two*QPrime (see Lehmer's definition).
void two_Qprime_in_range(long min, long max, char *res, mpz_t primes[], int num_primes);

int main(int argc, char **argv)
{
    // The number of Pell equations solved.
    long counter[NUM_THREADS];
    // Sieving time (searching for coefficients)
    clock_t sieving_time[NUM_THREADS];
    // Solving time (for solving pell equations)
    clock_t solving_time[NUM_THREADS];
    long start, end, step;
    printf("Starting coefficient: ");
    gmp_scanf("%ld", &start);
    printf("Ending coefficient: ");
    gmp_scanf("%ld", &end);
    printf("Step: ");
    gmp_scanf("%ld", &step);

    FILE *fp[NUM_THREADS];

    long cuts[NUM_THREADS + 1];
    long sub_interval = (end - start + (NUM_THREADS - 1)) / NUM_THREADS;
    cuts[0] = start;
    cuts[NUM_THREADS] = end;
    for (int i = 1; i < NUM_THREADS; i++)
    {
        cuts[i] = start + i * sub_interval;
    }

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    double wall_start_time = omp_get_wtime();

#pragma omp parallel num_threads(NUM_THREADS)
    {
        int thread = omp_get_thread_num();
        counter[thread] = 0;
        sieving_time[thread] = 0;
        solving_time[thread] = 0;

        char file_path[100];
        sprintf(file_path, "results/res_solve_pell_%d_%ld-%ld.txt", thread, cuts[thread], cuts[thread+1]);

        fp[thread] = fopen(file_path, "w");
        if (fp[thread] == NULL)
        {
            perror("Couldn't open file.");
            // return EXIT_FAILURE;
        }

        long num_steps = 0;
        long curr_start = cuts[thread];
        while (curr_start < cuts[thread + 1])
        {
            // printf(curr_start);
            sieve_interval_pell(curr_start, curr_start + step, fp[thread], b, primes, thread, &sieving_time[thread], &solving_time[thread], &counter[thread]);
            curr_start += step;
            num_steps++;
            if (num_steps % 1000 == 0)
            {
                printf("Steps complete (thread %d): %ld\n", thread, num_steps);
            }
        }

        printf("Equations solved (thread %d): %ld\n", thread, counter[thread]);

        fclose(fp[thread]);

        printf("Sieving Time (thread %d): %ld seconds\n", thread, sieving_time[thread]);
        printf("Solving Time (thread %d): %ld seconds\n", thread, solving_time[thread]);
    }

    printf("Total Time taken %lf seconds\n", omp_get_wtime() - wall_start_time);

    for (int i = 0; i < NUM_PRIMES; i++)
    {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;
}

void sieve_interval_pell(long start, long end, FILE *fp, mpz_t b, mpz_t primes[], int thread, clock_t *sieving_time, clock_t *solving_time, long *counter)
{
    clock_t sieve_start_time = 0;
    clock_t solve_start_time = 0;
    char res[end - start + 1]; // cover last pair to overlap
    for (long i = 0; i < end - start + 1; i++)
    {
        res[i] = 0;
    }
    sieve_start_time = omp_get_wtime();
    two_Qprime_in_range(start, end, res, primes, NUM_PRIMES);
    *sieving_time += omp_get_wtime() - sieve_start_time;
    for (long i = start; i < end; i++)
    {
        if (res[i - start] == 1)
        {
            (*counter)++;
            //printf("In 2QPRIME: %ld\n", (i));
            mpz_t result;
            mpz_init(result);
            mpz_t d;
            mpz_init(d);
            mpz_set_ui(d, i);
            if (mpz_cmp_ui(d, 4) != 0)
            { // not perfect square TODO Maybe replace with removing "4" from 2QPrime
                solve_start_time = omp_get_wtime();
                solve_pell(d, b, result, primes, NUM_PRIMES);
                *solving_time += omp_get_wtime() - solve_start_time;
                // gmp_printf("%Zd\n", result);
                if (mpz_cmp_si(result, 0) != 0)
                {
                    mpz_out_str(fp, 10, result);
                    fputs(" 1\n", fp);
                    check_and_compute_higher_solutions(result, primes, NUM_PRIMES, fp);
                }
            }
            mpz_clear(result);
            mpz_clear(d);
        }
    }
}

void two_Qprime_in_range(long min, long max, char *res, mpz_t primes[], int num_primes)
{ // including min, excluding max
    long nums[max - min];
    for (long i = 0; i < max - min; i++)
    {
        nums[i] = min + i;
    }
    for (int i = 0; i < num_primes; i++)
    {
        long div_by = mpz_get_ui(primes[i]);
        long mult = min / div_by + (min % div_by != 0);
        while (mult * div_by < max)
        {
            nums[mult * div_by - min] /= mpz_get_ui(primes[i]);
            mult += 1;
        }
    }
    for (long i = 0; i < max - min; i++)
    {
        // printf("start: %ld index: %ld value: %ld", min, i, nums[i]);
        if ((nums[i] == 2 || nums[i] == 1) && (i + min) % 2 == 0)
        {
            res[i] = 1;
        }
    }
}


