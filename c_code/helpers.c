#include <time.h>
#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

// number of primes below 32768
#define NUM_PRIMES 3512 
// smoothness bound
#define BOUND 32768
// The maximum period of continued fractions that we allow. 1000 should be very safe.
#define MAX_PERIOD 1000
// The cutoff for continued fractions. Will only find pairs up to this cutoff-1 bits from
// fundamental solutions.
#define BIT_CUTOFF (258)

// The number of Pell equations solved.
long counter;

// Sieving time (searching for coefficients)
clock_t sieving_time;

// Solving time (for solving pell equations)
clock_t solving_time;

void primes_up_to_b(mpz_t* ptr, mpz_t b) {
	int index = 0;
	mpz_t curr;
	mpz_init_set_str(curr, "2", 10);
	while (mpz_cmp(curr, b) < 0) {
		mpz_init_set(ptr[index], curr);
		mpz_nextprime(curr, curr);
		index++;
	}
	mpz_clear(curr);
}


bool is_smooth(mpz_t arg, mpz_t primes[], int num_primes) {
	mpz_t check;
	mpz_init(check);
	mpz_set(check, arg);
	for (int i = 0; i < num_primes; i++) {
		while (mpz_divisible_p(check, primes[i]) != 0) {
			mpz_divexact(check, check, primes[i]);
		}
	}
	mpz_t one;
	mpz_init_set_str(one, "1", 10);
	if (mpz_cmp(check, one) == 0) {
		mpz_clear(check);
		mpz_clear(one);
		return true;
	}
	mpz_clear(check);
	mpz_clear(one);
	return false;
}


void add_term(int coeff, int pow, mpz_t current, mpz_t eval) {
        mpz_t val;
        mpz_init(val);
        mpz_pow_ui(val, eval, pow);
        mpz_mul_si(val, val, coeff);
        mpz_add(current, val, current);
        mpz_clear(val);
}

// fourth polynomial
bool check_fourth_poly(mpz_t m, mpz_t primes[]) {
	mpz_t res1, res2;
	mpz_init(res1);
	mpz_add_ui(res1, res1, 1);
	add_term(2, 1, res1, m);
	//gmp_printf("First val: %Zd \n", res1);
	if (! is_smooth(res1, primes, NUM_PRIMES)) {
		mpz_clear(res1);
		return false;
	}
	mpz_init(res2);
	mpz_add_ui(res2, res2, 1);
	add_term(8, 1, res2, m);
	add_term(8, 2, res2, m);
	//gmp_printf("Second val: %Zd \n", res2);
	if (! is_smooth(res2, primes, NUM_PRIMES)) {
		mpz_clear(res1);
		mpz_clear(res2);
                return false;
        }
	mpz_clear(res1);
	mpz_clear(res2);
	return true;
}

bool check_second_poly(mpz_t m, mpz_t primes[]) {
    mpz_t res1;
    mpz_init(res1);
    mpz_add_ui(res1, res1, 1);
    add_term(2, 1, res1, m);
    if (! is_smooth(res1, primes, NUM_PRIMES)) {
	mpz_clear(res1);
        return false;
    }
    mpz_clear(res1);
    return true;
}

bool check_sixth_poly(mpz_t m, mpz_t primes[]) { // 2*(16*m^2 + 16*m + 1)*(4*m + 3)*(4*m + 1)*(2*m + 1)
    mpz_t res1, res2, res3, res4;
    mpz_init(res1);
    mpz_add_ui(res1, res1, 1);
    add_term(2, 1, res1, m);
    if (! is_smooth(res1, primes, NUM_PRIMES)) {
	mpz_clear(res1);
        return false;
    }
    mpz_init(res2);
    mpz_add_ui(res2, res2, 1);
    add_term(4, 1, res2, m);
    if (! is_smooth(res2, primes, NUM_PRIMES)) {
	mpz_clear(res1);
	mpz_clear(res2);
        return false;
    }
    mpz_init(res3);
    mpz_add_ui(res3, res3, 3);
    add_term(4, 1, res3, m);
    if (! is_smooth(res3, primes, NUM_PRIMES)) {
	mpz_clear(res1);
	mpz_clear(res2);
	mpz_clear(res3);
        return false;
    }
    mpz_init(res4);
    mpz_add_ui(res4, res4, 1);
    add_term(16, 1, res4, m);
    add_term(16, 2, res4, m);
    if (! is_smooth(res4, primes, NUM_PRIMES)) {
	mpz_clear(res1);
	mpz_clear(res2);
	mpz_clear(res3);
	mpz_clear(res4);
        return false;
    }
    mpz_clear(res1);
    mpz_clear(res2);
    mpz_clear(res3);
    mpz_clear(res4);
    return true;
}

bool check_seventh_poly(mpz_t m, mpz_t primes[]) { // (64*m^3 + 112*m^2 + 56*m + 7)*(64*m^3 + 80*m^2 + 24*m + 1)
    mpz_t res1, res2;
    mpz_init(res1);
    mpz_add_ui(res1, res1, 7);
    add_term(56, 1, res1, m);
    add_term(112, 2, res1, m);
    add_term(64, 3, res1, m);
    if (! is_smooth(res1, primes, NUM_PRIMES)) {
        mpz_clear(res1);
        return false;
    }
    mpz_init(res2);
    mpz_add_ui(res2, res2, 1);
    add_term(24, 1, res2, m);
    add_term(80, 2, res2, m);
    add_term(64, 3, res2, m);
    if (! is_smooth(res2, primes, NUM_PRIMES)) {
        mpz_clear(res1);
        mpz_clear(res2);
        return false;
    }
    mpz_clear(res1);
    mpz_clear(res2);
    return true;
}

void value_of_second_poly(mpz_t m, mpz_t result) {
    //mpz_add_ui(result, result, 1);
    add_term(4, 1, result, m);
    add_term(4, 2, result, m);
}

void value_of_fourth_poly(mpz_t m, mpz_t result) {  // 128*m^4 + 256*m^3 + 160*m^2 + 32*m + 1
    //mpz_add_ui(result, result, 1);
    add_term(16, 1, result, m);
    add_term(80, 2, result, m);
    add_term(128, 3, result, m);
    add_term(64, 4, result, m);
}


void smooths_in_range(mpz_t primes[], long min, long max, int num_primes, char* res) { // including min, excluding max
	long nums[max - min];
	for (long i = 0; i < max - min; i++) {
		nums[i] = min + i;
	}
	for (int i = 0; i < num_primes; i++) {
		long div_by = mpz_get_ui(primes[i]);
		while (div_by < max) {
			long mult = min/div_by + (min % div_by != 0);
			while (mult * div_by < max) {
				nums[mult * div_by - min] /= mpz_get_ui(primes[i]);
				mult += 1;
			}
		div_by *= mpz_get_ui(primes[i]);
		}
	}
	for (long i = 0; i < max - min; i++) {
		if (nums[i] == 1) {
			res[i] = 1;
		}
	}
}

bool check_pell(mpz_t x, mpz_t y, mpz_t d) {
    mpz_t xsquared;
    mpz_init(xsquared);
    mpz_mul(xsquared, x, x);
    mpz_t dysquared;
    mpz_init(dysquared);
    mpz_mul(dysquared, y, y);
    mpz_mul(dysquared, dysquared, d);
    mpz_sub(xsquared, xsquared, dysquared);
    if (mpz_cmp_ui(xsquared, 1) == 0) {
        mpz_clear(xsquared);
        mpz_clear(dysquared);
        return true;
    }
    mpz_clear(xsquared);
    mpz_clear(dysquared);
    return false;
}

void solve_pell(mpz_t d, mpz_t b, mpz_t result, mpz_t primes[], int num_primes) {
    // Uses the method of continued fractions, and the recurence relations on
    // page 382 of Rosen's book Elementary Number Theory to generate the convergents.
    // Cuts off when the numerator gets too high to save time.

    clock_t start_time = clock();
    
    mpz_t one;
    mpz_t zero;
    mpz_t cutoff;
    mpz_init_set_str(one, "1", 10);
    mpz_init_set_str(zero, "0", 10);
    mpz_init_set_str(cutoff, "1", 10);
    mpz_mul_2exp(cutoff, cutoff, BIT_CUTOFF);

    // lists:
    mpz_t numerators[MAX_PERIOD];
    mpz_t denominators[MAX_PERIOD];
    mpz_t p_k[MAX_PERIOD];
    mpz_t q_k[MAX_PERIOD];
    mpz_t a_k[MAX_PERIOD]; // convergents to the square root of d

    // intialization:
    mpz_init_set(p_k[0], zero);
    mpz_init_set(q_k[0], one);
    mpz_init(a_k[0]);
    mpz_sqrt(a_k[0], d);
    mpz_init_set(numerators[0], zero);
    mpz_init_set(numerators[1], one);
    mpz_init_set(numerators[2], a_k[0]);
    mpz_init_set(denominators[0], one);
    mpz_init_set(denominators[1], zero);
    mpz_init_set(denominators[2], one);
    int index = 2;

    mpz_t psquared;
    mpz_init(psquared);
    mpz_t p_plus_root_d;
    mpz_init(p_plus_root_d);

    // main loop
    while(!check_pell(numerators[index], denominators[index], d)) {
        if (mpz_cmp(numerators[index], cutoff) >= 0) {
            mpz_set(result, zero);
            for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
                mpz_clear(p_k[i]);
                mpz_clear(q_k[i]);
                mpz_clear(a_k[i]);
            }
            for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
                mpz_clear(numerators[i]);
                mpz_clear(denominators[i]);
            }
            mpz_clear(psquared);
            mpz_clear(p_plus_root_d);
	    mpz_clear(one);
            mpz_clear(zero);
            mpz_clear(cutoff);
            return; // this pell equation does not give useful smooths
        }
        mpz_init(p_k[index-1]);
        mpz_init(q_k[index-1]);
        mpz_init(a_k[index-1]);
        mpz_init(numerators[index+1]);
        mpz_init(denominators[index+1]);
        // generate the next convergent:
        // p_(k+1)
        // set_p_k_next(p_k[index-1], a_k[index-2], q_k[index-2], p_k[index-2]);
        mpz_mul(p_k[index-1], a_k[index-2], q_k[index-2]);
        mpz_sub(p_k[index-1], p_k[index-1], p_k[index-2]);
        // q_(k+1)
        // set_q_k_next(q_k[index-1], d, q_k[index-2], p_k[index-1]);
        mpz_mul(psquared, p_k[index-1], p_k[index-1]);
        mpz_sub(q_k[index-1], d, psquared);
        mpz_divexact(q_k[index-1], q_k[index-1], q_k[index-2]);
        // a_(k+1)
        // set_a_k_next(d, p_k[index-1], q_k[index-1]);
        // use the fact that floor(sqrt(d)) is a_k[0].
        mpz_add(p_plus_root_d, p_k[index-1], a_k[0]);
        mpz_fdiv_q(a_k[index-1], p_plus_root_d, q_k[index-1]);

        // generate new numerator and denominator:
        mpz_mul(numerators[index+1], numerators[index], a_k[index-1]);
        mpz_add(numerators[index+1], numerators[index+1], numerators[index-1]);

        mpz_mul(denominators[index+1], denominators[index], a_k[index-1]);
        mpz_add(denominators[index+1], denominators[index+1], denominators[index-1]);

        index++;
    }
    if (is_smooth(denominators[index], primes, num_primes)) {
        mpz_sub(result, numerators[index], one);
        mpz_divexact_ui(result, result, 2);
        solving_time += clock() - start_time;

        for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
            mpz_clear(p_k[i]);
            mpz_clear(q_k[i]);
            mpz_clear(a_k[i]);
        }
        for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
            mpz_clear(numerators[i]);
            mpz_clear(denominators[i]);
        }
        mpz_clear(psquared);
        mpz_clear(p_plus_root_d);
        mpz_clear(one);
        mpz_clear(zero);
        mpz_clear(cutoff);
        return;
    }
    mpz_set(result, zero);
    solving_time += clock() - start_time;

    for (int i = 0; i < index - 1; i++) { // clear ints in p_k, q_k, a_k
            mpz_clear(p_k[i]);
            mpz_clear(q_k[i]);
            mpz_clear(a_k[i]);
    }
    for (int i = 0; i < index + 1; i++) { // clear ints in numerators and denominators
        mpz_clear(numerators[i]);
        mpz_clear(denominators[i]);
    }
    mpz_clear(psquared);
    mpz_clear(p_plus_root_d);
    mpz_clear(one);
    mpz_clear(zero);
    mpz_clear(cutoff);
    return; // y was not smooth
}
