#include <time.h>
#include <omp.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "helpers.h"
#include "config.h"

int main(int argc, char **argv)
{
    // Sieving time
    double sieving_time[NUM_THREADS];
    // Time to check the higher solutions
    double checking_time[NUM_THREADS];

    long counter[NUM_THREADS];
    long start, end, step, min_bits;
    printf("Start: ");
    scanf("%ld", &start);
    printf("End: ");
    scanf("%ld", &end);
    printf("Step: ");
    scanf("%ld", &step);
    printf("Minimal output bitsize: ");
    scanf("%ld", &min_bits);

    // optimize finding twin smooths

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
        checking_time[thread] = 0;

        char file_path[200];
        sprintf(file_path, "results/res_sieve_and_check_%d_%ld-%ld.txt", thread, cuts[thread], cuts[thread + 1]);

        fp[thread] = fopen(file_path, "w");
        if (fp[thread] == NULL)
        {
            perror("Couldn't open file.");
            // return EXIT_FAILURE;
        }

        long curr_start = cuts[thread];
        while (curr_start < cuts[thread + 1])
        {
            char res[step + 1]; // cover last pair to overlap
            for (long i = 0; i < step + 1; i++)
            {
                res[i] = 0;
            }
            counter[thread] += 1;
            double sieve_start_time = omp_get_wtime();
            smooths_in_range(primes, curr_start, curr_start + step + 1, NUM_PRIMES, res);
            sieving_time[thread] += omp_get_wtime() - sieve_start_time;

            double check_start_time = omp_get_wtime();
            for (long i = 0; i < step; i++)
            {
                if ((res[i] == 1) && (res[i + 1] == 1))
                {
                    // Array of booleans for keeping track of whether 
                    // consecutive solutions correspond to smooth pairs.
                    bool is_pair[13];                  
                    // The first number m in the smooth pair (m, m+1).
                    mpz_t m;
                    mpz_init_set_si(m, i + curr_start);
                    // Bitsize of m.
                    size_t m_bits = mpz_sizeinbase(m, 2);
                    // printf("SMOOTH: %ld \n ", i + curr_start);
                    
                    // Check wether higher consecutive solutions also give pairs.
                    check_higher_solutions(m, primes, NUM_PRIMES, is_pair);
                    int n = 1;
                    while ((n*m_bits < min_bits) && (n < 13)) {
                        n++;
                    }
                    while (!is_pair[n] && (n < 13)) {
                        n++;
                    }
                    if (n < 13) {
                        mpz_out_str(fp[thread], 10, m);
                        for (int j = n; j < 13; j++) {
                            if (is_pair[j]) {
                                fprintf(fp[thread], " %d", j);
                            }
                        }
                        fputs("\n", fp[thread]);
                    }
                
                    mpz_clear(m);
                }
            }
            checking_time[thread] += omp_get_wtime() - check_start_time;

            curr_start += step;

            if (counter[thread] % 100 == 0)
            {
                printf("[%d] -> %ld, %ld intervals, sieving: %lfs, checking: %lfs\n",
                       thread, curr_start + step, counter[thread], sieving_time[thread], checking_time[thread]);
            }
        }

        printf("Total sieving Time (thread %d): %lf seconds\n", thread, sieving_time[thread]);
        printf("Total checking Time (thread %d): %lf seconds\n", thread, checking_time[thread]);

        fclose(fp[thread]);
    }

    printf("Total Time taken %lf seconds\n", omp_get_wtime() - wall_start_time);

    for (int i = 0; i < NUM_PRIMES; i++)
    {
        mpz_clear(primes[i]);
    }
    mpz_clear(b);

    return EXIT_SUCCESS;
}
