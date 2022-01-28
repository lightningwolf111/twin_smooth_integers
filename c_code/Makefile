PROGS = solve_pell factor_method pell_sols_4_6_7 read_check_2_4 read_max solve_pell_parallel connected_component_search pell_sols_4_6_7_parallel small_primes_search small_primes_search_parallel

all :$(PROGS)

solve_pell: solve_pell.c helpers.c helpers.h
	gcc -g -o solve_pell solve_pell.c helpers.c -lgmp

solve_pell_parallel: solve_pell_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o solve_pell_parallel solve_pell_parallel.c helpers.c -lgmp

small_primes_search: small_primes_search.c helpers.c helpers.h
	gcc -g -o small_primes_search small_primes_search.c helpers.c -lgmp

small_primes_search_parallel: small_primes_search_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o small_primes_search_parallel small_primes_search_parallel.c helpers.c -lgmp

factor_method: factor_method.c helpers.c helpers.h
	gcc -g -o factor_method factor_method.c helpers.c -lgmp

pell_sols_4_6_7: pell_sols_4_6_7.c helpers.c helpers.h
	gcc -g -o pell_sols_4_6_7 pell_sols_4_6_7.c helpers.c -lgmp

pell_sols_4_6_7_parallel: pell_sols_4_6_7_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o pell_sols_4_6_7_parallel pell_sols_4_6_7_parallel.c helpers.c -lgmp

read_check_2_4: read_check_2_4.c helpers.c helpers.h
	gcc -g -o read_check_2_4 read_check_2_4.c helpers.c -lgmp

read_max: read_max.c helpers.c helpers.h
	gcc -g -o read_max read_max.c helpers.c -lgmp

connected_component_search: connected_component_search.c helpers.c helpers.h solve_pell_extended.cc solve_pell_extended.h 
	g++ -g -o connected_component_search connected_component_search.c helpers.c solve_pell_extended.cc -lgmp -lgmpxx

conrey_near_pell: conrey_near_pell.cc helpers.c helpers.h solve_pell_extended.cc solve_pell_extended.h
	g++ -g -o conrey_near_pell conrey_near_pell.cc solve_pell_extended.cc helpers.c -lgmp -lgmpxx

clean:
	rm -f *.o $(PROGS) fourth_pell_solution
