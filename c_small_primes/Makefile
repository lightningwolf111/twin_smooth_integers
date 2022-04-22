PROGS = read_max small_primes_search_parallel

all :$(PROGS)

small_primes_search_parallel: small_primes_search_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o small_primes_search_parallel small_primes_search_parallel.c helpers.c -lgmp

read_max: read_max.c helpers.c helpers.h
	gcc -g -o read_max read_max.c helpers.c -lgmp

clean:
	rm -f *.o $(PROGS)
