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

void readFrom(char* read_from, char* write_2, char* write_4, mpz_t* primes, mpz_t b);

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

    if (S_ISDIR(root_stat.st_mode)) {
        DIR *rd;
        rd = opendir(filename);
        struct dirent *dirent;
	char write_2[100];
	char write_4[100];
        for (int i = 0; (dirent = readdir(rd)) != NULL; i++) {
            if (strcmp(dirent->d_name, ".") != 0 &&
                strcmp(dirent->d_name, "..") != 0) {
                char *to_read;
                strcpy(to_read, filename);
                to_read = strcat(to_read, dirent->d_name);
                //printf("Reading from %s\n", to_read);
		//char write_2[100];
                // don't include .txt
                strncpy(write_2, to_read, strlen(to_read) - 4);
                strcat(write_2, "_2nd_sol.txt");
                //char write_4[100];
                strncpy(write_4, to_read, strlen(to_read) - 4);
                strcat(write_4, "_4th_sol.txt");
                readFrom(to_read, write_2, write_4, primes, b);
		printf("filename: %s\n", filename);
		printf("to_read: %s\n", to_read);
		printf("write_2: %s\n", write_2);
		printf("write_4 %s\n", write_4);
		memset(write_2, '\0', sizeof write_2); // resets temp
		memset(write_4, '\0', sizeof write_4); // resets temp
		memset(to_read, '\0', sizeof to_read);
            }
        }
    } else {  // we assume it is a file
        char write_2[100];
        // don't include .txt
        strncpy(write_2, filename, strlen(filename) - 4);
        strcat(write_2, "_2nd_sol.txt");
        char write_4[100];
        strncpy(write_4, filename, strlen(filename) - 4);
        strcat(write_4, "_4th_sol.txt");
        readFrom(filename, write_2, write_4, primes, b);
    }


    mpz_clear(b);
    for (int i = 0; i < NUM_PRIMES; i++) {
            mpz_clear(primes[i]);
    }
    return EXIT_SUCCESS;
}


void readFrom(char* filename, char* write_2, char* write_4, mpz_t* primes, mpz_t b) {
    FILE *read_from;
    FILE *write_to_2;
    FILE *write_to_4;

    write_to_2 = fopen(write_2, "w");
    write_to_4 = fopen(write_4, "w");
    read_from = fopen(filename, "r");

    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    mpz_t m;
    mpz_init(m);
    mpz_t second;
    mpz_t fourth;
    while ((read = getline(&line, &len, read_from)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        mpz_set_str(m, line, 10);
        //gmp_printf("m: %Zd\n", m);
        if (check_second_poly(m, primes)) {
            mpz_init(second);
            value_of_second_poly(m, second);
            mpz_out_str(write_to_2, 10, second);
            fputs("\n", write_to_2);
            mpz_clear(second);
        }
        if (check_fourth_poly(m, primes)) {
            mpz_init(fourth);
            value_of_fourth_poly(m, fourth);
            mpz_out_str(write_to_4, 10, fourth);
            fputs("\n", write_to_4);
            mpz_clear(fourth);
        }
    }

    mpz_clear(m);

    fclose(read_from);
    fclose(write_to_2);
    fclose(write_to_4);
}
