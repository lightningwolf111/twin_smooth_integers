#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include "helpers.h"

// number of primes below 32768
#define NUM_PRIMES 3512
// smoothness bound
#define BOUND 32768

void readFrom(char* filename, mpz_t* primes, mpz_t b, mpz_t* max_ret);

int main(int argc, char **argv) {
    char filename[100];
    gmp_printf("Filename: \n");
    scanf("%s", filename);

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    struct stat root_stat;
    if (stat((char *) filename, &root_stat) == -1) {
        printf("Error retrieving file stats");
    }

    mpz_t max;
    mpz_t max_current;
    mpz_init(max);
    mpz_init(max_current);

    if (S_ISDIR(root_stat.st_mode)) {
        DIR *rd;
        rd = opendir(filename);
        struct dirent *dirent;
        for (int i = 0; (dirent = readdir(rd)) != NULL; i++) {
	if (strcmp(dirent->d_name, ".") != 0 &&
            strcmp(dirent->d_name, "..") != 0) {
		char to_read[100];
		strcpy(to_read, filename);
		strcat(to_read, dirent->d_name);
		//printf("Reading from %s\n", to_read);
                readFrom(to_read, primes, b, &max_current);
		if (mpz_cmp(max, max_current) < 0) {
                    mpz_set(max, max_current);
		}
		printf("%d\n", i);
            }
        }
	gmp_printf("Overall max: %Zd\n", max);
    } else {  // we assume it is a file
        readFrom(filename, primes, b, &max);
    }
    mpz_clear(max);
    mpz_clear(max_current);
    mpz_clear(b);
    for (int i = 0; i < NUM_PRIMES; i++) {
            mpz_clear(primes[i]);
    }
    return EXIT_SUCCESS;
}


void readFrom(char* filename, mpz_t* primes, mpz_t b, mpz_t* max_ret) {
    FILE *read_from;

    read_from = fopen(filename, "r");   
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    mpz_t m;
    mpz_init(m);
    mpz_t max;
    mpz_init(max);
    while ((read = getline(&line, &len, read_from)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        mpz_set_str(m, line, 10);
        if (mpz_cmp(m, max) > 0) {
            mpz_set(max, m);
        }
    }
    if (mpz_sizeinbase(max, 2) > 180) {
    	gmp_printf("Maximum pair: %Zd\n in file %s\n", max, filename);
    }
    mpz_set(*max_ret, max);

    mpz_clear(m);
    mpz_clear(max);

    fclose(read_from);
}
