#Adjusted Conrey code to generate b(B+1) or b(B+1)B(b+1) observations (square free parts)
#Assumes the existence of four folders: Full_b(B+1)_data/ and Full_b(B+1)B(b+1)_data/ to store raw data, and
#Full_b(B+1)_observations/ and Full_b(B+1)B(b+1)_observations/ to store the observations.

import math
from itertools import combinations

def primesUpToB(b):
    primes = []
    P = Primes()
    primes.append(P.first())
    next = P.next(2)
    while next <= b:
        primes.append(next)
        next = P.next(next)
    return primes

def square_free(i):
    ret = []
    factors = factor(i)
    for f in factors:
        if f[1] % 2 == 1:
           ret.append(f[0])
    string = str(ret).replace(",", "")
    return string[1:len(string) - 1]

def find_beta(b, B, file, func):
    fraction = (b*(B+1)/((b+1)*B))
    if (fraction.numerator() + 1) == fraction.denominator():
        file.write(square_free(func(b, B)) + " \n")
        return fraction.numerator()
    return -1

def inQ(primes, y):
    for prime in primes:
        while (y % prime == 0):
            y = y / prime
    if (y == 1):
        return True
    else:
        return False
    
def expand_set_one_iteration_opt(S, new_elements, file, func):
    to_add = set()
    to_check = set(combinations(new_elements, 2))
    for element in new_elements:
        for s in S:
                if s < element:
                        to_check.add((s, element))
                elif s > element:
                        to_check.add((element, s))

    for element in to_check:
        res = find_beta(element[0], element[1], file, func)
        if (not (res == -1)) and (not res in S): #removed "and (not res in new_elements)" condition
            to_add.add(res)

    return (to_add.union(S), to_add)

def expand_set_to_max_opt(S, file, func):
    result = expand_set_one_iteration_opt(S, S, file, func)
    s_prime = result[0]
    added = result[1]
    while not len(added) == 0:
        S = s_prime
        result = expand_set_one_iteration_opt(S, added, file, func)
        s_prime = result[0]
        added = result[1]
    return s_prime


# Optimized verion that only checks for new pairs.
# Runs Conrey's method starting with the numbers up to B.
def conrey_most_smooths_opt(b, num_fact):
    S = set([i + 1 for i in range(b)])
    if (num_fact == 2):
        file = open("Full_b(B+1)_data/" + str(b) + "_smooth_b(B+1)_data.txt", 'x')
        func = lambda b, B: b * (B+1)
    else:
        file = open("Full_b(B+1)B(b+1)_data/" + str(b) + "_smooth_b(B+1)B(b+1)_data.txt", 'x')
        func = lambda b, B: b * B * (b+1) * (B+1)
    result = expand_set_to_max_opt(S, file, func)
    res_list = list(result)
    res_list.sort()
    return res_list

def generate_factorizations(primelist, dict):
    for j in range(len(primelist)):
        for i in range(len(primelist) - 1):
            for tuple in combinations(primelist[j:len(primelist)], i + 1):
                dict[str(list(tuple))] = 0
    dict[str(primelist)] = 0


#run with num_fact = 2 to get b(B+1) observations, or num_fact = 4 to get b(B+1)B(b+1) observations
def generate_bB_observations_experimental(prime_list, num_fact):
    if 2 in prime_list:
        prime_list.remove(2)
    for prime in prime_list:
        conrey_most_smooths_opt(prime, num_fact)
        if (num_fact == 2):
            file = open("Full_b(B+1)_data/" + str(prime) + "_smooth_b(B+1)_data.txt", 'r')
        else:
            file = open("Full_b(B+1)B(b+1)_data/" + str(prime) + "_smooth_b(B+1)B(b+1)_data.txt", 'r')
        lines = file.readlines()
        primes = primesUpToB(prime)
        
        dict = {}
        generate_factorizations(primes, dict)
        
        for line in lines:
            arr = line.split()
            for i in range(len(arr)):
                dict[str([int(x) for x in arr[0:i + 1]])] += 1
            
        if (num_fact == 2):
            write_file = open("Full_b(B+1)_observations/" + str(prime) + "_smooth_b(B+1)_observations.txt", 'x')
        else:
            write_file = open("Full_b(B+1)B(b+1)_observations/" + str(prime) + "_smooth_b(B+1)B(b+1)_observations.txt", 'x')
            
        for key in dict.keys():
            write_file.write(str(key) + ": " + str(dict[key]) + "\n")
