PROGS = solve_pell solve_pell_parallel factor_method pell_sols_4_6_7

all :$(PROGS)

solve_pell: solve_pell.c helpers.c helpers.h
	gcc -g -o solve_pell solve_pell.c helpers.c -lgmp

solve_pell_parallel: solve_pell_parallel.c helpers.c helpers.h
	gcc -g -fopenmp -o solve_pell_parallel solve_pell_parallel.c helpers.c -lgmp

factor_method: factor_method.c helpers.c helpers.h
	gcc -g -o factor_method factor_method.c helpers.c -lgmp

pell_sols_4_6_7: pell_sols_4_6_7.c helpers.c helpers.h
	gcc -g -o pell_sols_4_6_7 pell_sols_4_6_7.c helpers.c -lgmp

clean:
	rm -f *.o $(PROGS) fourth_pell_solution
