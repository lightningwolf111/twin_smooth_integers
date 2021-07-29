#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

// number of primes below 32000
#define NUM_PRIMES 3432 
// smoothness bound
#define BOUND 32000

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
