def read_fundamental_solutions(filename, cutoff, sum_prime_required, check_fourth_sols):
    file = open(filename, 'r')
    L = [l for l in file.readlines()]
    nums = [((Integer) (l)) for l in L]
    nums.sort()
    file.close()
    print(len(nums))
    r = generate_up_to_3(7)
    primes = primesUpToB(32000)
    print("")
    print("")
    print("FUNDAMENTAL SOLUTIONS")
    print("")
    print("")
    for sm in nums:
        if (sm > 2**cutoff and (not sum_prime_required or is_prime(sm*2+1))):
            print(sm)
            print(factor(sm))
            print(factor(sm+1))
            print(ceil(log(sm, 2)))
            print(is_prime(sm*2+1))
    print("")
    print("")
    print("SECOND SOLUTIONS")
    print("")
    print("")
    for sm in nums:
        if (sm > 2**(cutoff/2) and inQ(primes, (2*sm)+1)):
            result_2 = r[0][2].subs(m = ((Integer) (sm)))
            smooth = ((Integer) ((result_2-1)/2))
            print(smooth)
            print(factor(smooth))
            print(factor(smooth+1))
            print(ceil(log(smooth, 2)))
            print(is_prime(result_2))
    print("")
    print("")
    print("FOURTH SOLUTIONS")
    print("")
    print("")
    if (check_fourth_sols):
        for sm in nums:
            if (sm > 2**(cutoff/4) and inQ(primes, (8*sm*sm + 8*sm + 1)) and inQ(primes, (2*sm+1))):
                result_4 = r[0][4].subs(m = ((Integer) (sm)))
                smooth = ((Integer) ((result_4-1)/2))
                print(smooth)
                print(factor(smooth))
                print(factor(smooth+1))
                print(ceil(log(smooth, 2)))
                print(is_prime(result_4))
