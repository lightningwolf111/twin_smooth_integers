# Returns a solution for the negative Pell equation with coefficient D.
def solve_neg_pell(D):
    if (D == 1):
        return -1
    reg_sol = fundamentalPellSolConvergentsExperimental(D)
    if (reg_sol == -1):
        return -1
    x_val = 0
    x_vals = solve(2*x^2 + 1 == reg_sol[0], x, solution_dict=True)
    for val in x_vals:
        x_val = max(val[x], x_val)
    if (not x_val in ZZ):
        return -1
    y_val = reg_sol[1]/(2*x_val)
    if (not y_val in ZZ):
        return -1
    return (x_val, y_val)


# Given a list of primes, returns those that are either 1 mod 4 or equal 2.
def primes_one_mod_four(primes):
    res = []
    for prime in primes:
        if (prime % 4 == 1 or prime == 2):
            res.append(prime)
    return res

# Returns a list of solutions to the negative pell equations with coefficients
# that are b-smooth and squarefree.
def solve_neg_pell_equations(b):
    sols = []
    reduced_primes = primes_one_mod_four(primesUpToB(b))
    for D in generateQPrime(reduced_primes):
        sol = solve_neg_pell(D)
        if (not sol == -1):
            sols.append((D, sol[0], sol[1]))
    return sols

# Returns the m from smooth pairs that result from solving the negative pell equations.
def smooths_from_neg(b):
    primes = primesUpToB(b)
    pairs = []
    sols = solve_neg_pell_equations(b)
    for sol in sols:
        x = ((Integer) (sol[1]))
        if inQ(primes, x) and inQ(primes, x**2 + 1):
            pairs.append(x**2)
    pairs.sort()
    return pairs
