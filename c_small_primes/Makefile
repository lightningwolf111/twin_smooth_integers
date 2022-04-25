PROGS = read_max solve_pell_parallel small_primes_search_parallel

all :$(PROGS)

solve_pell_parallel: solve_pell_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o solve_pell_parallel solve_pell_parallel.c helpers.c -lgmp

small_primes_search_parallel: small_primes_search_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o small_primes_search_parallel small_primes_search_parallel.c helpers.c -lgmp

read_max: read_max.c helpers.c helpers.h
	gcc -g -o read_max read_max.c helpers.c -lgmp

clean:
	rm -f *.o $(PROGS)
