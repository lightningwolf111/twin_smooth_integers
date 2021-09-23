'''
THIS CODE DOES NOT WORK.
'''

from sage.modules.free_module_integer import IntegerLattice

def primesUpToB(b):
    primes = []
    P = Primes()
    primes.append(P.first())
    next = P.next(2)
    while next <= b:
        primes.append(next)
        next = P.next(next)
    return primes

def gamma(d):
  return sqrt(d) * (2**(d/2))

def a_opt(d, B, P):
  return d * sqrt ( (((1 - (1/d)) * log(B, 2))/(log(P, 2))) + ((2 * log( gamma(d), 2))/(d * log(P, 2))) )

def epsilon(d, B, P):
  return (sqrt( ( (log(B, 2))/ (log(P, 2)) ) + ( 1 - (1/d) ) + ( ((2 * log(gamma(d), 2))) / (d * log(P, 2)) )) + ( 5 / (4 * d)))

def g(i, a, P, R, B):
  t = ((P**(a-i)) * (((x*B) - R)**i))
  return t.coefficients(sparse=False)

def h(i, a, R, B):
  t = (((x*B) - R) ** a) * ((x*B) ** i)
  return t.coefficients(sparse=False)

def I_bound(T, S, d):
  expon = ( (log(T, 2) / log(S, 2)) - (2.5/d) )
  return (T ** expon)/4

def amp(m, P, R):
  return gcd(P, m-R)

def crt_decode(B, p_list, r_list, d):
  # required constants:
  P = 1
  for p in p_list:
    P *= p

  # this will NOT work for general CRT decoding. It's a shortcut since we know all r_i = -U
  R = P + r_list[0]

  a = round(a_opt(d, B, P))
  a_prime = d - a

  # step 1: construct a polynomial w(x) over the integers so that all required m are roots of w(x)

  # step 1a: given d > 3 construct a dxd matrix whose rows are the coeffficients of the g and h polynomials in xB,
  # with a being the closest integer to a_opt
  L_arr = []

  # generate the g family of polynomials
  for i in range(0, a):
    coeffs_short = g(i, a, P, R, B) #without trailing zeroes
    coeffs_full = coeffs_short + [0 for i in range(0, d - len(coeffs_short))] #with trailing zeroes
    L_arr.append(coeffs_full)

  # generate the h family of polynomials
  for i in range(0, a_prime):
    coeffs_short = h(i, a, R, B) #without trailing zeroes
    coeffs_full = coeffs_short + [0 for i in range(0, d - len(coeffs_short))] #with trailing zeroes
    L_arr.append(coeffs_full)

  # step 1b: Run the LLL algorithm on the d-dimensional lattice spanned by the
  # rows of the matrix L. Let v be the resulting short vector in L.

  print("B=%s\nR=%s\nP=%s" %(B, R, P))
  #print(matrix(L_arr))
  L = IntegerLattice(L_arr, lll_reduce=True)
  v = L.shortest_vector()

  # step 1c: View v as the coefficients of a d−1 degree polynomial w(xB).
  # Output w(x) as the required polynomial
  #v = v[::-1]
  w(x) = 0
  w_coeffs = [val / (B**i) for i, val in enumerate(v)]
  for i, val in enumerate(w_coeffs):
    w += val * (x**(i))

  #print(w)
  # step 2: find all integer roots of w(x)
  amp_bound = P ** (sqrt( ( log (4 * B, 2)/ log(P, 2)) ) + (5 / (4*d)))

  rts = [k[0] for k in w.roots() if (k[0] in ZZ and amp(k[0], P, R) > amp_bound)]

  solns = [r - r_list[0] for r in rts]

  print("Found solutions:")
  print(solns)

  return solns

def crt_strongly_smooth(s, T, d, abs_I):
  '''
  uses crt decoding to find all strongly s-smooth
  integers in the interval [2T - abs_I, 2T]
  '''

  pr = primesUpToB(s+1)
  S = 1
  P = 1

  q_list = []
  for q in pr:
    a_i =  floor(log(s, 2) / log(q, 2))
    S *= q ** (a_i)
    q_list += [q ** a_i]
    P *= q

  min_d = 1 + sqrt ( log(S, 2) / log(abs_I, 2) )
  max_width = I_bound(T, S, d)
  if d < min_d:
    print("Lattice dimension too small. Minimum with (s = %s, interval width = 2**%s) is %s." %(s, log(abs_I, 2).n(20), min_d.n(20)))
  elif abs_I > max_width:
    print("Interval too large. Maximum with (s = %s, T = 2**%s, d = %s) is ~ 2**%s." %(s, log(T, 2).n(20), d, log(max_width, 2).n(20)))
  else:
    U = 2 * T - floor(abs_I)
    V = 2 * T
    u_list = [-U for i in range(0, len(q_list))]

    # instance of crt problem is now given by abs_I, q_list, and u_list
    return crt_decode(abs_I, q_list, u_list, d)