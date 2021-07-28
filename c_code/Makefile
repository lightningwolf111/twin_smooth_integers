all: solve_pell factor_method fourth_pell_solution

solve_pell: solve_pell.c helpers.c helpers.h
	gcc -g -o solve_pell solve_pell.c helpers.c -lgmp 

factor_method: factor_method.c helpers.c helpers.h
	gcc -g -o factor_method factor_method.c helpers.c -lgmp

fourth_pell_solution: fourth_pell_solution.c helpers.c helpers.h
	gcc -g -o fourth_pell_solution fourth_pell_solution.c helpers.c -lgmp

clean:
	rm -f *.o solve_pell factor_method fourth_pell_solution
