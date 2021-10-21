import time


def fundamentalPellSolConvergentsUnoptimized(d):
    cont_frac = continued_fraction(sqrt(d))
    frac = (cont_frac[0] + 1 / cont_frac[1])
    index = 2
    while not (frac.numerator()**2 - frac.denominator()**2 * d == 1):
        # print(str(frac.numerator()**2 - frac.denominator()**2 * d == 1))
        currIndex = index
        currF = cont_frac[index]
        while currIndex > 0:
            currIndex = currIndex - 1
            currF = cont_frac[currIndex] + 1 / currF
        # print(currF)
        frac = currF
        index = index + 1
    return (frac.numerator(), frac.denominator())

def write_Pell_fund_to_file(n):
    sols = []
    for i in range(n):
        if i % 100 == 0:
            print("Solved " + str(i) + " equations")
        if not Integer(2*i).is_square():
            sols.append(fundamentalPellSolConvergentsUnoptimized(2*i))
        else:
            sols.append((0, 0))
    file = open("pell_sols.txt", 'w')
    file.write(str(sols))
    file.close()

def write_Pell_dictionary_to_file(b):
    sols = {}
    QPrime = generateQPrime(primesUpToB(b))
    QPrime.remove(2)
    for num in QPrime:
        sols[num*2] = fundamentalPellSolConvergentsUnoptimized(num * 2)
    file = open("pell_sols_dict.txt", 'w')
    file.write(str(sols))
    file.close()

def read_Pell_fund_from_file():
    file = open("pell_sols.txt", 'r')
    s = file.readline()
    file.close()
    sols = []
    sols.extend(eval(s))
    return sols

def read_Pell_dictionary_from_file():
    file = open("pell_sols_dict.txt", 'r')
    dictionary = file.readline()
    file.close()
    return eval(dictionary)




pell_solutions = []
pell_dictionary = {}

# load("Write_Pell_Solutions_File.sage")
print("loading sols")
pell_solutions = read_Pell_fund_from_file()
print("loading dictionary")
# load("Write_Pell_Solutions_File.sage")
pell_dictionary = read_Pell_dictionary_from_file()


# copied from bounds:
def primesUpToB(b):
    primes = []
    P = Primes()
    primes.append(P.first())
    next = P.next(2)
    while next <= b:
        primes.append(next)
        next = P.next(next)
    return primes

def findM(primes):
    return max(3, (primes[len(primes)-1] + 1)/2)

def inQ(primes, y):
    for prime in primes:
        while (y % prime == 0):
            y = y / prime
    if (y == 1):
        return True
    else:
        return False

def fundamentalPellSol(d):
    x = 2
    while not sqrt((x**2 - 1) / d) in ZZ:
        x = x + 1
    return (x, sqrt((x**2 - 1) / d))

# A copy of the above, taking advantage of the fact that in our case, x in solutions is odd.
def fundamentalPellSolOpt(d):
    x = 3
    while not sqrt((x**2 - 1) / d) in ZZ:
        x = x + 2
    return (x, sqrt((x**2 - 1) / d))

def fundamentalPellSolConvergents(d):
    global pell_solutions
    global pell_dictionary
    if pell_solutions == []:
        # load("Write_Pell_Solutions_File.sage")
        print("loading sols")
        pell_solutions = read_Pell_fund_from_file()
    if d % 2 == 0 and d < len(pell_solutions*2):
        # print("Optimized!")
        return pell_solutions[d/2]
    
    if pell_dictionary == {}:
        print("loading dictionary")
        # load("Write_Pell_Solutions_File.sage")
        pell_dictionary = read_Pell_dictionary_from_file()
    if d in pell_dictionary.keys():
        # print("dictionary used!")
        return pell_dictionary[d]
    
    print("Not found")
    cont_frac = continued_fraction(sqrt(d))
    frac = (cont_frac[0] + 1 / cont_frac[1])
    index = 2
    while not (frac.numerator()**2 - frac.denominator()**2 * d == 1):
        # print(str(frac.numerator()**2 - frac.denominator()**2 * d == 1))
        currIndex = index
        currF = cont_frac[index]
        while currIndex > 0:
            currIndex = currIndex - 1
            currF = cont_frac[currIndex] + 1 / currF
        # print(currF)
        frac = currF
        index = index + 1
    return (frac.numerator(), frac.denominator()) 

def nextPellSol(x, y, xfund, yfund, d):
    nextSol =  ((x + sqrt(d)*y) * (xfund + yfund*sqrt(d))).full_simplify()
    tupleForm = (nextSol.operands())
    return (tupleForm[1], tupleForm[0] / sqrt(d))

def generateQPrime(primes):
    QPrime = []
    if len(primes) == 0:
        return [1]
    add = primes.pop()
    for num in generateQPrime(primes):
        QPrime.append(num*add)
        QPrime.append(num)
    return QPrime

def solvePellEquations(QPrime, M, b):
    solutions = []
    counter = 1
    for num in QPrime:
        print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
        fundSol = fundamentalPellSolConvergents(2*num)
        solutions.append(fundSol)
        newSol = fundSol
        # print("Generating extra solutions")
        if inQ(primesUpToB(b), Integer(newSol[1])):
            for index in range(2, M+1):
                newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
                solutions.append(newSol)
        counter = counter + 1;
    return solutions

def findSmooth(b):
    primes = primesUpToB(b)
    M = findM(primes)
    QPrime = generateQPrime(primes)
    QPrime.remove(2)
    # print("QPrime: " + str(QPrime))
    solutions = solvePellEquations(QPrime, M, b)
    smoothPairs = []
    for sol in solutions:
        if inQ(primesUpToB(b), Integer(sol[1])):
            smoothNum = (sol[0] - 1) / 2
            # print("Solution " + str(smoothNum))
            smoothPairs.append(smoothNum)
    smoothPairs.sort()
    print(str(smoothPairs[-1]))

# finds smooths and information about which Pell solutions generated pairs.
# returns a list of 3-tuples where
# first element is smoothNum, next is coeff. of equation divided by 2, next is which num solution it is.
def findSmoothWithInfo(b):
    primes = primesUpToB(b)
    M = findM(primes)
    QPrime = generateQPrime(primes)
    QPrime.remove(2)
    sols = [] # first element is smoothNum, next is coeff. of equation divided by 2, next is which num solution it is.
    counter = 1
    for num in QPrime:
        print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
        fundSol = fundamentalPellSolConvergents(2*num)
        newSol = fundSol
        if inQ(primesUpToB(b), Integer(fundSol[1])):
            smoothNum = (fundSol[0] - 1) / 2
            sols.append((smoothNum, num, 1))
            # print("Fundamental Solution gives" + str(smoothNum))
        for index in range(2, M+1):
            newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
            if inQ(primesUpToB(b), Integer(newSol[1])):
                smoothNum = (newSol[0] - 1) / 2
                sols.append((smoothNum, num, index))
                # print(str(index) + "th power Solution gives " + str(smoothNum))
        counter = counter + 1;
    sols.sort()
    return sols


def solvePellEquationsFaster(QPrime, M, b):
    solutions = []
    counter = 1
    for num in QPrime:
        print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
        fundSol = fundamentalPellSolConvergents(2*num)
        solutions.append(fundSol)
        newSol = fundSol
        print("Generating extra solutions")
        if inQ(primesUpToB(b), Integer(newSol[1])):
            for index in range(2, 30):
                newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
                solutions.append(newSol)
        counter = counter + 1;
    return solutions



def fundamentalPellSolConvergentsExperimental(d):
    st = time.time()
    global pell_solutions
    global pell_dictionary
    if pell_solutions == []:
        # load("Write_Pell_Solutions_File.sage")
        # print("loading sols")
        pell_solutions = read_Pell_fund_from_file()
    if d % 2 == 0 and d < len(pell_solutions*2):
        # print("Optimized!")
        return pell_solutions[d/2]
    
    if pell_dictionary == {}:
        print("loading dictionary")
        # load("Write_Pell_Solutions_File.sage")
        pell_dictionary = read_Pell_dictionary_from_file()
    if d in pell_dictionary.keys():
        # print("dictionary used!")
        return pell_dictionary[d]
    
    #print("Not found")
    cont_frac = continued_fraction(sqrt(d))
    frac = (cont_frac[0] + 1 / cont_frac[1])
    index = 2
    while not (frac.numerator()**2 - frac.denominator()**2 * d == 1):
        # print(str(frac.numerator()**2 - frac.denominator()**2 * d == 1))
        currIndex = index
        currF = cont_frac[index]
        while currIndex > 0:
            currIndex = currIndex - 1
            currF = cont_frac[currIndex] + 1 / currF
        # print(currF)
        frac = currF
        index = index + 1
        if (frac.numerator() > 2**260):
            return -1
    # print(time.time() - st)
    return (frac.numerator(), frac.denominator())


def generateQPrimeExperimental(primes, max):
    QPrime = []
    if len(primes) == 0:
        return [1]
    add = primes.pop()
    for num in generateQPrimeExperimental(primes, max):
        if (num*add < max):
            QPrime.append(num*add)
        QPrime.append(num)
    return QPrime



def findSmoothWithInfoExperimental(b, max):
    primes = primesUpToB(b)
    extra_primes = primesUpToB(b)
    M = findM(primes)
    QPrime = generateQPrimeExperimental(primes, max)
    QPrime.remove(2)
    sols = [] # first element is smoothNum, next is coeff. of equation divided by 2, next is which num solution it is.
    # print("Number of equations: " + str(len(QPrime)))
    counter = 1
    for num in QPrime:
        #print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
        fundSol = fundamentalPellSolConvergentsExperimental(2*num)
        if (not (fundSol == -1)):
            newSol = fundSol
            if inQ(extra_primes, Integer(fundSol[1])):
                smoothNum = (fundSol[0] - 1) / 2
                sols.append((smoothNum, num, 1))
         #       print("Fundamental Solution gives" + str(smoothNum))
            #for index in range(2, M+1):
            #    newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
            #    if inQ(primesUpToB(b), Integer(newSol[1])):
            #        smoothNum = (newSol[0] - 1) / 2
            #        sols.append((smoothNum, num, index))
            #        # print(str(index) + "th power Solution gives " + str(smoothNum))
        counter = counter + 1;
    sols.sort()
    return sols


import numpy as np

def randomExperimental(b, num_coeff_facts, quantity, skew, max_b):
    primes = primesUpToB(b)
    extra_primes = primesUpToB(max_b)
    counter = 1
    distribution = []
    for i in range(len(primes)):
        distribution.append((1/primes[i])**skew)
    sum = 0
    for i in distribution:
        sum += i
    for i in range(len(primes)):
        distribution[i] = distribution[i]/sum
    # print(distribution)
    factors = np.random.multinomial(num_coeff_facts, distribution, size=quantity)
    QPrime = set()
    sols = [] # first element is smoothNum, next is coeff. of equation divided by 2, next is which num solution it is.
    for arr in factors:
        coeff = 1
        for i in range(len(arr)):
            coeff *= (primes[i])^(arr[i])
        if (not sqrt(2*coeff) in ZZ):
            QPrime.add(coeff)
        # print(coeff)
    for num in QPrime:
        # print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
        fundSol = fundamentalPellSolConvergentsExperimental(2*num)
        if (not (fundSol == -1 or fundSol[0] == 0)):
            newSol = fundSol
            if inQ(extra_primes, Integer(fundSol[1])):
                smoothNum = (fundSol[0] - 1) / 2
                sols.append((smoothNum, num, 1))
                print("Fundamental Solution gives " + str(smoothNum))
            #for index in range(2, M+1):
            #    newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
            #    if inQ(primesUpToB(b), Integer(newSol[1])):
            #        smoothNum = (newSol[0] - 1) / 2
            #        sols.append((smoothNum, num, index))
            #        # print(str(index) + "th power Solution gives " + str(smoothNum))
        counter = counter + 1;
        if counter % 100 == 0:
            print(counter)
    sols.sort()
    return sols

def solvePellFully(coeff):
    extra_primes = primesUpToB(2**15)
    M = findM(extra_primes)
    fundSol = fundamentalPellSolConvergentsExperimental(2*coeff)
    if (not (fundSol == -1 or fundSol[1] == 0)):
        newSol = fundSol
        print(fundSol)
        if inQ(extra_primes, Integer(fundSol[1])):
            smoothNum = (fundSol[0] - 1) / 2
            # sols.append((smoothNum, coeff, 1))
            print("Fundamental Solution gives " + str(smoothNum))
            for index in range(2, M+1):
                print(index) 
                newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*coeff)
                if inQ(extra_primes, Integer(newSol[1])):
                    smoothNum = (newSol[0] - 1) / 2
                    # sols.append((smoothNum, coeff, index))
                    print(str(index) + "th power Solution gives " + str(smoothNum))

def solveMostly(coeff, b):
    extra_primes = primesUpToB(b)
    M = findM(extra_primes)
    fundSol = fundamentalPellSolConvergentsExperimental(2*coeff)
    if (not (fundSol == -1 or fundSol[1] == 0)):
        newSol = fundSol
        print(fundSol)
        if inQ(extra_primes, Integer(fundSol[1])):
            smoothNum = (fundSol[0] - 1) / 2
            # sols.append((smoothNum, coeff, 1))
            print("Fundamental Solution gives " + str(smoothNum))
            for index in range(2, 12):
                # print(index)
                newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*coeff)
                if inQ(extra_primes, Integer(newSol[1])):
                    smoothNum = (newSol[0] - 1) / 2
                    # sols.append((smoothNum, coeff, index))
                    print(str(index) + "th power Solution gives " + str(smoothNum))

def gives_smooth(coeff, b):
    primes = primesUpToB(b)    
    fundSol = fundamentalPellSolConvergentsExperimental(2*squarefree_part(coeff))
    if (not (fundSol == -1 or fundSol[1] == 0)):
        if inQ(primes, Integer(fundSol[1])):
            return true
    return false

def pell_testing(b):
    primes = primesUpToB(b)
    num_not = 0
    num_not_pq = 0
    num_not_for_p_yes = 0
    for p in primes:
        if (not gives_smooth(p, b)):
            num_not += 1
            for q in primes:
                if (not gives_smooth(p*q, b)):
                    num_not_pq += 1
        else:
            for q in primes:
                if (not gives_smooth(p*q, b)):
                    num_not_for_p_yes += 1
    print("Yes for p: " + str(len(primes) - num_not))
    print("No for pq when yes for p: " + str(num_not_for_p_yes))
    print("No for p: " + str(num_not))
    print("No for pq when no for p: " + str(num_not_pq))
    print("Ratio for p yes: " + str(1 - (num_not_for_p_yes /(len(primes) * (len(primes) - num_not)))))
    print("Ratio for p no: " + str(1 - (num_not_pq /(len(primes) * (num_not)))))


def fundamentalPellSolConvergentsExperimental_2(d):
    # st = time.time()
    cont_frac = continued_fraction(sqrt(d))
    numerators = [0, 1, cont_frac[0]]
    denominators = [1, 0, 1]
    index = 2
    while not (numerators[index]**2 - denominators[index]**2 * d == 1):
        index = index + 1
        numerators.append(cont_frac[index - 2] * numerators[index - 1] + numerators[index - 2])
        denominators.append(cont_frac[index - 2] * denominators[index - 1] + denominators[index - 2])
        print("Numerator at " + str(index) + " is " + str(numerators[index]))
        print("Denominator at " + str(index) + " is " + str(denominators[index]))
        print("Convergent: " + str(cont_frac[index - 2]))
        if (numerators[index] > 2**500):
            return -1
    # print(time.time() - st)
#    print(index)
    return (numerators[index], denominators[index])

def pellWhichConvergent(d):
    d = squarefree_part(d)
    if d == 1:
        return -1
    cont_frac = continued_fraction(sqrt(d))
    numerators = [0, 1, cont_frac[0]]
    denominators = [1, 0, 1]
    index = 2
    while not (numerators[index]**2 - denominators[index]**2 * d == 1):
        index = index + 1
        numerators.append(cont_frac[index - 2] * numerators[index - 1] + numerators[index - 2])
        denominators.append(cont_frac[index - 2] * denominators[index - 1] + denominators[index - 2])
        if (numerators[index] > 2**300):
            return -1
    return index

def whichConvergentsInRange(min, max):
    unsolved = 0
    solutions = [0 for i in range(200)]
    for i in smooths_in_range(32000, min, max):
        res = pellWhichConvergent(i + min)
        if res == -1:
            unsolved += 1
        else:
            solutions[res] += 1
    print("Unsolved: " + str(unsolved))
    return solutions


def average(l):
    sum = 0
    nums = 0
    for i in range(len(l)):
        sum += i*l[i]
        nums += l[i]
    return (sum / nums).numerical_approx()

def numInRange(min, max, dist):
    i = 2 ** min
    while (i < 2**max):
        num_solved = 0
        num_work = 0
        smooths = smooths_in_range(32000, i, i + 2**dist)
        for smooth in smooths:
            if (squarefree_part(smooth) == smooth):
                num_solved += 1
                res = fundamentalPellSolConvergentsExperimental(smooth)
                if (res != -1 and res[0] > 240):
                    num_work += 1
        print("Range starts at: " + str(i) + " num solved: " + str(num_solved) + " over 240 bits " + str(num_work) + " fraction " + str((num_work / num_solved).numerical_approx()))
        i *= 2




def findSmoothWithInfoExperimental_2(b, min, max):
    primes = primesUpToB(b)
    extra_primes = primesUpToB(32000)
    M = findM(primes)
    QPrime = generateQPrimeExperimental(primes, max)
    QPrime.remove(2)
    sols = [] # first element is smoothNum, next is coeff. of equation divided by 2, next is which num solution it is.
    print("Number of equations: " + str(len(QPrime)))
    counter = 1
    for num in QPrime:
        if num > min:
            print("SolvingPellEq number " + str(counter) + " with coefficient of " + str(2*num))
            fundSol = fundamentalPellSolConvergentsExperimental(2*num)
            if (not (fundSol == -1)):
                newSol = fundSol
                if inQ(extra_primes, Integer(fundSol[1])):
                    smoothNum = (fundSol[0] - 1) / 2
                    sols.append((smoothNum, num, 1))
                    print("Fundamental Solution gives" + str(smoothNum))
            #for index in range(2, M+1):
            #    newSol = nextPellSol(newSol[0], newSol[1], fundSol[0], fundSol[1], 2*num)
            #    if inQ(primesUpToB(b), Integer(newSol[1])):
            #        smoothNum = (newSol[0] - 1) / 2
            #        sols.append((smoothNum, num, index))
            #        # print(str(index) + "th power Solution gives " + str(smoothNum))
        counter = counter + 1;
    sols.sort()
    return sols

def which_frac(d):
    # st = time.time()
    cont_frac = continued_fraction(sqrt(d))
    numerators = [0, 1, cont_frac[0]]
    denominators = [1, 0, 1]
    index = 2
    while not (numerators[index]**2 - denominators[index]**2 * d == 1):
        index = index + 1
        numerators.append(cont_frac[index - 2] * numerators[index - 1] + numerators[index - 2])
        denominators.append(cont_frac[index - 2] * denominators[index - 1] + denominators[index - 2])
        # print("Numerator at " + str(index) + " is " + str(numerators[index]))
        # print("Denominator at " + str(index) + " is " + str(denominators[index]))
        if (numerators[index] > 2**260):
            return -1
    # print(time.time() - st)
    print(index)
    return (index)

# Returns a list of tuples, in which the first element is prime <= b,
# and the second is the number of times that it occurs in coefficients
# of pell equations that do not give any smooth pairs.
def pell_not_working(b):
    res = [0 for i in range(b+1)]
    primes = primesUpToB(b)
    QPrime = generateQPrime(primes)
    QPrime.remove(2)
    primes = primesUpToB(b)
    for coeff in QPrime:
        sol = fundamentalPellSolConvergentsExperimental(2*coeff);
        if ((sol == -1) or (not inQ(primes, sol[1]))):
            for fact in factor(coeff):
                res[fact[0]] += 1
    return res

from itertools import  combinations

# Returns a list of tuples, where the first two elements are primes,
# and the last is how many times these lead to coefficients with no
# smooths resulting from the pell equation.
def pell_bad_pairs(b):
    primes = primesUpToB(b)
    res = {}
    for comb in combinations(primes, 2):
        res[(comb[0], comb[1])] = 0
    QPrime = generateQPrime(primes)
    QPrime.remove(2)
    primes = primesUpToB(b)
    for coeff in QPrime:
        sol = fundamentalPellSolConvergentsExperimental(2*coeff);
        if ((sol == -1) or (not inQ(primes, sol[1]))):
            # print("coeff: " + str(coeff) + " y: " + str(sol[1]))
            factor_set = [i[0] for i in factor(coeff)]
            for pair in combinations(factor_set, 2):
                res[pair] += 1
    return res

# Generalization of the above to n-tuples
def pell_bad_n_tuples(n, b):
    primes = primesUpToB(b)
    res = {}
    for comb in combinations(primes, n):
        res[comb] = 0
    QPrime = generateQPrime(primes)
    QPrime.remove(2)
    primes = primesUpToB(b)
    for coeff in QPrime:
        sol = fundamentalPellSolConvergentsExperimental(2*coeff);
        if ((sol == -1) or (not inQ(primes, sol[1]))):
            # print("coeff: " + str(coeff) + " y: " + str(sol[1]))
            factor_set = [i[0] for i in factor(coeff)]
            for pair in combinations(factor_set, n):
                res[pair] += 1
    return res

###########################################################################

# Investigating certain particular classes of pell equations.

def get_factors(current, p_facts):
    if (len(p_facts) == 0):
        return [current]
    res = []
    prime_tup = p_facts[0]
    for i in range(prime_tup[1] + 1):
        newList = list(p_facts)
        newList.pop(0)
        res += get_factors(current * prime_tup[0] ** i, newList)
    return res


def generate_factors(Q):
    fact = factor(4*Q)
    return get_factors(1, list(fact))



# Solve Pell equations of Richaud-Degert type for positive R:
def solve_pell_RD_type_pos(minQ, maxQ):
    maxSol = 0
    minSol = 2**1000000
    primes = primesUpToB(2**15)
    smooths = []
    for Q in range(minQ, maxQ):
        for R in generate_factors(Q):
            D = Q**2 + R
#            D = squarefree_part(D)
            if (D % 2 == 1):
                D = D * 4
            if (D > 1 and D % 2 == 0):
                sol = fundamentalPellSolConvergentsExperimental_2(D)
                if (sol == -1):
                    print("NOT FOUND")
                    continue;
                if (sol[0] > maxSol):
                    maxSol = sol[0]
                print("Q: " + str(Q) + " R: " + str(R) + " D: " +  str(D) + "   Solution: " + str(sol))
                if inQ(primes, Integer(sol[1])):
                    smoothNum = (sol[0] - 1) / 2
                    smooths.append(smoothNum)
    print("Maximum solution: " + str(maxSol))
    return smooths



# Solve Pell equations of Richaud-Degert type for negative R:
def solve_pell_RD_type_neg(minQ, maxQ):
    maxSol = 0
    minSol = 2**1000000
    smooths = []
    primes = primesUpToB(2**15)
    for Q in range(minQ, maxQ):
        for R in generate_factors(Q):
            D = Q**2 - R
            D = squarefree_part(D)
            if (D % 2 == 1):
                D = D * 4
            if (D > 1 and D % 2 == 0):
                sol = fundamentalPellSolConvergentsExperimental_2(D)
                if (sol == -1):
                    print("NOT FOUND")
                    continue;
                if (sol[0] > maxSol):
                    maxSol = sol[0]
                print("Q: " + str(Q) + " R: " + str(R) + " Solution: " + str(sol))
                if inQ(primes, Integer(sol[1])):
                    smoothNum = (sol[0] - 1) / 2
                    smooths.append(smoothNum)
    print("Maximum solution: " + str(maxSol))
    return smooths


# Experimental, solves the above for only R = 1.
def solve_pell_RD_type_pos_1(minQ, maxQ):
    maxSol = 0
    minSol = 2**1000000
    primes = primesUpToB(2**15)
    smooths = []
    for Q in range(minQ, maxQ):
        for R in [1]:
            D = Q**2 + R
#            D = squarefree_part(D)
            if (D % 2 == 1):
                D = D * 4
            if (D > 1 and D % 2 == 0):
                sol = fundamentalPellSolConvergentsExperimental_2(D)
                if (sol == -1):
                    print("NOT FOUND")
                    continue;
                if (sol[0] > maxSol):
                    maxSol = sol[0]
                print("Q: " + str(Q) + " R: " + str(R) + " D: " +  str(D) + "   Solution: " + str(sol))
                if inQ(primes, Integer(sol[1])):
                    smoothNum = (sol[0] - 1) / 2
                    smooths.append(smoothNum)
    print("Maximum solution: " + str(maxSol))
    return smooths


def generate_factors_without_times_4(Q):
    fact = factor(Q)
    return get_factors(1, list(fact))

# An algorithm based on D = Q^2 + R type pell equations that sieves for Q and looks for suitable R.
def Q_and_R_algorithm(minQ, maxQ, b):
    smooths = []
    primes = primesUpToB(b)
    for Q in range(minQ, maxQ):
        if (not inQ(primes, Q)):
            continue;
        R_values = generate_factors_without_times_4(Q)
        R_values.sort()
        for R in R_values:
            if (inQ(primes, Q**2 + R)):
                smooths.append(Q**2/R)
                continue;
    return smooths

# Returns any smooth numbers in the given range.
def find_smooth_in_range(min, max, b):
    smooths = []
    primes = primesUpToB(b)
    for i in range(min, max):
        if (inQ(primes, i)):
            smooths.append(i)
    return smooths

