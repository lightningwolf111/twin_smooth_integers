import numpy as np
import time

# The start time
st = 0

# Generates a 'random' coefficient between min and min*b.
def generate_random(min, primes):
    num = 2
    while (num < min):
        p = primes[floor(random() * len(primes))]
        if (p == 2):
            if (num % 4 == 0):
                continue
        else:
            if (num % p == 0):
                continue
        num *= p
    return num

# Generates a random squarefree smooth between min and min * b, according to the given distribution.
def generate_random_skewed(min, primes, distribution):
    ints = [i for i in range(len(primes))]
    distribution = list(distribution) # Make a local copy for modification
    num = 2
    primesChosen = np.random.choice(ints, 50, False, distribution)
    index = 0
    while (num < min):
        p = primes[primesChosen[index]]
        if num % p != 0:
             num *= p
        index += 1
        if index == 50:
            primesChosen = np.random.choice(ints, 50, False, distribution)
            index = 0
    return num

# Like the above, but instead begins with the number being 'base'.
def generate_random_with_base(base, min, primes, distribution):
    global st
    method_st = time.time()
    ints = [i for i in range(len(primes))]
    distribution = list(distribution) # Make a local copy for modification
    num = base
    primesChosen = np.random.choice(ints, 20, False, distribution)
    index = 0
    while (num < min):
        p = primes[primesChosen[index]]
        if num % p != 0:
             num *= p
        index += 1
        if index == 20:
            primesChosen = np.random.choice(ints, 20, False, distribution)
            index = 0
    st += time.time() - method_st
    return num

def generate_skewed_distribution(primes, skew):
    distribution = []
    for i in range(len(primes)):
        distribution.append((1/primes[i])**skew)
    sum = 0
    for i in distribution:
        sum += (i).numerical_approx()
    for i in range(len(primes)):
        distribution[i] = (distribution[i]/sum).numerical_approx()
    return distribution


def generate_log_skewed_distribution(primes, skew):
    distribution = []
    for i in range(len(primes)):
        distribution.append((1/log(primes[i]))**skew)
    sum = 0
    for i in distribution:
        sum += i.numerical_approx()
    for i in range(len(primes)):
        distribution[i] = (distribution[i]/sum).numerical_approx()
    return distribution


def compare_distributions(min, dist1, dist2, trials, primes):
    sols1 = []
    sols2 = []
    for i in range(trials):
        num1 = generate_random_skewed(min, primes, dist1)
        num2 = generate_random_skewed(min, primes, dist2)
        res1 = fundamentalPellSolConvergentsExperimental(2*num1)
        res2 = fundamentalPellSolConvergentsExperimental(2*num2)
        if (res1 != -1):
            sols1.append(res1)
        if (res2 != -1):
            sols2.append(res2)
        if (i % 100 == 0):
            print("Finished " + str(i))
    sols1 = set(sols1)
    sols2 = set(sols2)
    print("Chance of getting a result for distribution 1: " + str(len(sols1)/trials))
    print("Chance of getting a result for distribution 2: " + str(len(sols2)/trials))
    return (sols1, sols2)

# Compares distributions, with the second having the given base.
def compare_distributions_base(min, dist1, dist2, trials, primes, base):
    sols1 = []
    sols2 = []
    for i in range(trials):
        num2 = generate_random_with_base(base, min, primes, dist1)
        num1 = generate_random_skewed(min, primes, dist2)
        res1 = fundamentalPellSolConvergentsExperimental(2*num1)
        res2 = fundamentalPellSolConvergentsExperimental(2*num2)
        if (res1 != -1):
            sols1.append(res1)
        if (res2 != -1):
            sols2.append(res2)
        if (i % 100 == 0):
            print("Finished " + str(i))
    sols1 = set(sols1)
    sols2 = set(sols2)
    print("Chance of getting a result for distribution 1: " + str(len(sols1)/trials))
    print("Chance of getting a result for distribution 2: " + str(len(sols2)/trials))
    return (sols1, sols2)

# As the above, but just with one distribution
def eval_distribution(min, dist1, trials, primes, base):
    global st
    sols = []
    for i in range(trials):
        num = generate_random_with_base(base, min, primes, dist1)
        res = fundamentalPellSolConvergentsExperimental_2(2*num)
        if (res != -1):
            sols.append(res)
        if (i % 100 == 0):
            print("Finished " + str(i))
    print("Number of results for distribution with repeats: " + str(len(sols)))
    sols = set(sols)
    print("Number of results for distribution: " + str(len(sols)))
    print("Time: " + str(st))
    st = 0
    return sols

# Disallow repeatted coefficients.
def eval_distribution_no_repeats(min, dist1, trials, primes, base):
    coeffs = set()
    global st
    sols = []
    for i in range(trials):
        num = generate_random_with_base(base, min, primes, dist1)
        while num in coeffs:
            num = generate_random_with_base(base, min, primes, dist1)
        res = fundamentalPellSolConvergentsExperimental_2(2*num)
        if (res != -1):
            sols.append(res)
        if (i % 100 == 0):
            print("Finished " + str(i))
        coeffs.add(num)
    print("Number of results for distribution with repeats: " + str(len(sols)))
    sols = set(sols)
    print("Number of results for distribution: " + str(len(sols)))
    print("Time: " + str(st))
    st = 0
    return sols



# Computes the statistical z-score for the hypothesis that p1 = p2. Largers absolute values
# mean that the true ps are 'more likely' different. 
def z_score(p1, p2, trials):
    p = (p1 + p2)/2
    return ((p1 - p2) / sqrt(p * (1-p) * (2/trials))).numerical_approx()

# rescales the given distribution so that values sum to 1.
def rescale(dist):
    sum = 0
    for i in dist:
        sum += (i).numerical_approx()
    for i in range(len(dist)):
        dist[i] = (dist[i]/sum).numerical_approx()
    return dist

###############################################################################################################
# Comparison with conrey


# num_generators is the number of working coefficients to generate the new conrey set from.
# gen_min is the minimum size of a generator (should be >= sqrt(min))
def eval_conrey(min, trials, primes, dist, num_generators, gen_min):
    base = 2
    coeffs = set()
    global st
    sols = []
    coeffs_half_bits_plus = [] # The conrey coefficients from which to generate the new ones
    while len(coeffs_half_bits_plus) < num_generators:
        num = generate_random_with_base(base, gen_min, primes, dist)
        while num in coeffs_half_bits_plus:
            num = generate_random_with_base(base, gen_min, primes, dist)
        res = fundamentalPellSolConvergentsExperimental_2(2*num)
        if (res != -1 and inQ(primes, res[1])):
            coeffs_half_bits_plus.append(num)
#            print("Building from " + str(num))
    print("Finished generating half size coefficients")
    for i in range(trials):
        num = 0
        while num < min:
            first = randrange(num_generators)
            second = randrange(num_generators)
            num = coeffs_half_bits_plus[first] * coeffs_half_bits_plus[second]
            num = squarefree_part(num)
#        print("Example: " + str(num))
        res = fundamentalPellSolConvergentsExperimental_2(2*num)
        if (res != -1): 
            sols.append(res)


    print("Number of results for distribution with repeats: " + str(len(sols)))
    sols = set(sols)
    print("Number of results for distribution: " + str(len(sols)))
    print("Time: " + str(st))
    st = 0
    return sols













###############################################################################################################
# Counting squarefree smooths (what is the distribution of primes in D?)

# Returns the index of the prime in the list, or -1 if not found
def whichIndex(primes, num):
    for i in range(len(primes)):
        if num == primes[i]:
            return i
    return -1


# Returns the squarefree smooths in the given interval
def squarefree_smooths_in_range(min, max, primes):
    res = []
    for i in range(max - min):
        if (squarefree_part(i + min) == i + min) and inQ(primes, i + min):
            res.append(i + min)
        if (i % 10000 == 0):
            print(i)
    return res


# Takes a list of squarefree numbers and returns the normalized vector of the primes in them:
def prime_frequencies(primes, squarefrees):
    dist = [0 for i in range(len(primes))]
    for num in squarefrees:
        fact = factor(num)
        for i in range(len(fact)):
            dist[whichIndex(primes, fact[i][0])] += 1
    st = time.time()
    return rescale(dist)



###############################################################################################################
# Comparing distributions: Wasserstein distance

def wasserstein_distance(dist1, dist2):
    if not len(dist1) == len(dist2):
        print("Distributions must have equal length")
        return
    sum = 0
    currentDiff = 0
    for i in range(len(dist1)):
        currentDiff += dist1[i] - dist2[i]
        sum += currentDiff
    return (sum/len(dist1)).numerical_approx()
