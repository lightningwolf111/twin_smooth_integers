# File for trying to get more precise bounds on the number of squarefree coefficients D
# that are b-smooth and less than K. Currently, the best obvious estimate is the Dickman
# function, but that does not restrict us to squarefree numbers.

# Helper function for generate_ys, which only uses values from start_index onward.
# current is the current vector y being created.
# log_k measures the max sum from start_index onward
def generate_ys_from(delta, counts, log_k, start_index, current):
    #print(start_index)
    #print(current)
    #print(log_k)
    ys = []
    if (start_index >= len(counts)):  # finished recursing
        ys.append(current)
        return ys
    for num in range(0, counts[start_index] + 1):  #include right endpoint
        if (log_k - num*delta*(start_index + 1) >= 0):  # we would still be under log_k
            to_send = list(current)
            to_send.append(num)
            res_list = generate_ys_from(delta, counts, log_k - num*delta*(start_index + 1), start_index+1, list(to_send))
            for res in res_list:
                ys.append(res)
    return list(ys)



# Generates all vectors y of the same dimension as counts such that the sum over i
# of i * delta * y_i is less than or equal to log(k).
def generate_ys(delta, counts, log_k):    
    return generate_ys_from(delta, counts, log_k, 0, [])


# Returns the sum over ys of the product of (count[i] choose y[i]).
# with generate_ys, this can bound the number of coefficients D.
def how_many_sums(ys, counts):
    sum = 0;
    for y in ys:
        prod = 1
        for i in range(0, len(counts)):
            prod *= binomial(counts[i], y[i])
        sum += prod
    return sum

# Does the exact calculation of the number of squarefree D less than 2**log_k.
# Slow/intractable for large inputs; useful for testing
def exact_squarefree_coeffs(b, log_k):
    count = 0
    for i in range(1, 2**log_k + 1):
        if (i == 1):
            count += 1
            continue
        fact = factor(i)
        if (fact[-1][0] <= b):
            toAdd = True
            for j in range(len(fact)):
                # print(str(fact) + " j: " + str(j))
                if (fact[j][1] != 1):
                    toAdd = false
            if (toAdd):
                #print(i)
                count += 1;
    return count
