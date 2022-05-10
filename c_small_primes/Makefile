PROGS = read_max solve_pell_parallel small_primes_search_parallel sieve_and_check_parallel

all :$(PROGS)

sieve_and_check_parallel: sieve_and_check_parallel.c helpers.c helpers.h config.h
	gcc -g -fopenmp -o sieve_and_check_parallel sieve_and_check_parallel.c helpers.c -lgmp

solve_pell_parallel: solve_pell_parallel.c helpers.c helpers.h config.h
	gcc -g -fopenmp -o solve_pell_parallel solve_pell_parallel.c helpers.c -lgmp

small_primes_search_parallel: small_primes_search_parallel.c helpers.c helpers.h config.h
	gcc -g -fopenmp -o small_primes_search_parallel small_primes_search_parallel.c helpers.c -lgmp

read_max: read_max.c helpers.c helpers.h
	gcc -g -o read_max read_max.c helpers.c -lgmp

clean:
	rm -f *.o $(PROGS)
