#include <time.h>
#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "helpers.h"

// The number of solutions in the desired range, regardless of whether y was smooth or not
long numInRange;

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

static void add_term(mpz_t current, int coeff, mpz_t eval, int pow) {
    // Computes current += coeff*eval^pow
    mpz_t val;
    mpz_init(val);
    mpz_pow_ui(val, eval, pow);  // val = eval^pow
    mpz_mul_si(val, val, coeff); // val = val*coeff = coeff*eval^pow
    mpz_add(current, val, current); // current += coeff*eval^pow
    mpz_clear(val);
}

// fourth polynomial
bool check_fourth_poly(mpz_t m, mpz_t primes[]) {
	mpz_t res1, res2;
	mpz_init(res1);
	mpz_add_ui(res1, res1, 1);
	add_term(res1, 2, m, 1); // res += 2*m^1
	//gmp_printf("First val: %Zd \n", res1);
	if (! is_smooth(res1, primes, NUM_PRIMES)) {
		mpz_clear(res1);
		return false;
	}
	mpz_init(res2);
	mpz_add_ui(res2, res2, 1);
	add_term(res2, 8, m, 1);
	add_term(res2, 8, m, 2);
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
    add_term(res1, 2, m, 1);
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
    add_term(res1, 2, m, 1);
    if (! is_smooth(res1, primes, NUM_PRIMES)) {
	mpz_clear(res1);
        return false;
    }
    mpz_init(res2);
    mpz_add_ui(res2, res2, 1);
    add_term(res2, 4, m, 1);
    if (! is_smooth(res2, primes, NUM_PRIMES)) {
	mpz_clear(res1);
	mpz_clear(res2);
        return false;
    }
    mpz_init(res3);
    mpz_add_ui(res3, res3, 3);
    add_term(res3, 4, m, 1);
    if (! is_smooth(res3, primes, NUM_PRIMES)) {
	mpz_clear(res1);
	mpz_clear(res2);
	mpz_clear(res3);
        return false;
    }
    mpz_init(res4);
    mpz_add_ui(res4, res4, 1);
    add_term(res4, 16, m, 1);
    add_term(res4, 16, m, 2);
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
    add_term(res1, 56, m, 1);
    add_term(res1, 112, m, 2);
    add_term(res1, 64, m, 3);
    if (! is_smooth(res1, primes, NUM_PRIMES)) {
        mpz_clear(res1);
        return false;
    }
    mpz_init(res2);
    mpz_add_ui(res2, res2, 1);
    add_term(res2, 24, m, 1);
    add_term(res2, 80, m, 2);
    add_term(res2, 64, m, 3);
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
    add_term(result, 4, m, 1);
    add_term(result, 4, m, 2);
}

void value_of_fourth_poly(mpz_t m, mpz_t result) {  // 128*m^4 + 256*m^3 + 160*m^2 + 32*m + 1
    //mpz_add_ui(result, result, 1);
    add_term(result, 16, m, 1);
    add_term(result, 80, m, 2);
    add_term(result, 128, m, 3);
    add_term(result, 64, m, 4);
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
    // Uses the method of continued fractions, and the recurrence relations on
    // page 382 of Rosen's book Elementary Number Theory to generate the convergents.
    // Cuts off when the numerator gets too high to save time.

    clock_t start_time = clock();
    clock_t solving_time = 0;
    
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
	    solving_time += clock() - start_time;
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

    ////////////////////
    mpz_t minBound;
    mpz_init_set_str(minBound, "1", 10);
    mpz_mul_2exp(minBound, minBound, BIT_CUTOFF - 18); // Set 2^240 bits as the min
    if (mpz_cmp(numerators[index], minBound) > 0) {
        numInRange++;
    }
    mpz_clear(minBound);
    ////////////////////

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


void square_free_part(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes) {
    mpz_set_ui(result, 1);
    mpz_t check;
    mpz_init(check);
    mpz_set(check, arg);
    for (int i = 0; i < num_primes; i++) {
        bool odd = false;
        while (mpz_divisible_p(check, primes[i]) != 0) {
            odd = !odd;
            mpz_divexact(check, check, primes[i]);
        }
	if (odd) {
            mpz_mul(result, result, primes[i]);
        }
    }

}

void extract_D(mpz_t arg, mpz_t result, mpz_t primes[], int num_primes) {
    mpz_t res_m_plus_one;
    mpz_init(res_m_plus_one);
    mpz_t m_plus_one;
    mpz_init(m_plus_one);
    mpz_add_ui(m_plus_one, arg, 1);
    square_free_part(arg, result, primes, num_primes);
    square_free_part(m_plus_one, res_m_plus_one, primes, num_primes);
    mpz_mul(result, result, res_m_plus_one);
    if (mpz_divisible_ui_p(arg, 4) != 0) {
        mpz_mul_ui(result, result, 4);
    }
}


void generate_file_with_coefficients(char* file_name) {
    FILE* read_from = fopen(file_name, "r");
    char write_to_string[100];
    strncpy(write_to_string, file_name, strlen(file_name) - 4);
    strcat(write_to_string, "_with_coeffs.txt");
    FILE* write_to = fopen(write_to_string, "w");

    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    mpz_t m;
    mpz_init(m);
    mpz_t d;
    mpz_init(d);

    mpz_t primes[NUM_PRIMES];
    mpz_t b;
    mpz_init_set_si(b, BOUND);
    primes_up_to_b(primes, b);

    while ((read = getline(&line, &len, read_from)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        mpz_set_str(m, line, 10);
        //gmp_printf("m: %Zd\n", m);
	extract_D(m, d, primes, NUM_PRIMES);
        mpz_out_str(write_to, 10, d);
	fputs(" ", write_to);
	mpz_out_str(write_to, 10, m);
        fputs("\n", write_to);
    }


    mpz_clear(m);
    mpz_clear(d);

    fclose(read_from);
    fclose(write_to);
}



void check_and_compute_higher_solutions(mpz_t m, mpz_t primes[], int num_primes, FILE *outputfile) {
    // Array of booleans for keeping track of whether 
    // consecutive solutions correspond to smooth pairs.
    bool is_pair[13];
    // The first two are not necessary, but included to have:
    // is_pair[i] = true if i-th solutions corresponds
    // to a smooth pair, else is_pair[i] = false.
    is_pair[0] = true;   
    is_pair[1] = true;
    for (int i=2; i<11; i++) {
        is_pair[i] = false;
    }
    // Initialize space for m2 as it might be used for
    // computing others such as m6, m10 and m12, not purely in the 
    // power-of-2 chain.
    mpz_t m2;
    mpz_init(m2);
    // Initialize space for wi0, wi1 for i in {3,5,7,11}. They are always 
    // computed to check for smoothness of pairs corresponding to the i-th
    // solutions.
    mpz_t w30, w31, w50, w51, w70, w71, w110, w111;
    mpz_init(w30);
    mpz_init(w31);
    mpz_init(w50);
    mpz_init(w51);
    mpz_init(w70);
    mpz_init(w71);
    mpz_init(w110);
    mpz_init(w111);
    
    // Compute x = 2*m+1 and its square x^2.
    mpz_t x, x2;
    mpz_init(x);
    mpz_init(x2);
    mpz_set_si(x, 1);
    mpz_addmul_ui(x, m, 2);
    mpz_mul(x2, x, x);
    
    // Now check smoothness of higher solutions by computing the wi.
    
    // Check 2nd solution:
    // corresponds to smooth pair if w2 = x is smooth.
    if (is_smooth(x, primes, num_primes)){
        mpz_t w4;
        mpz_init(w4);

        is_pair[2] = true;
        // Compute the pair (m2, m2+1) corresponding to the 2nd solution.
        mpz_add_ui(m2, m, 1);   // m2 = m+1
        mpz_mul(m2, m2, m);     // m2 = m*(m+1)
        mpz_mul_ui(m2, m2, 4);  // m2 = 4*m*(m+1)
        mpz_out_str(outputfile, 10, m2);
        fputs(" 2", outputfile);
        fputs("\n", outputfile);
        
        // Check 4th solution:
        // corresponds to smooth pair if w4 = 2*x^2-1 is smooth.
        mpz_add(w4, x2, x2);
        mpz_sub_ui(w4, w4, 1);
        if (is_smooth(w4, primes, num_primes)){
            mpz_t m4, w8;
            mpz_init(w8);
            mpz_init(m4);

            is_pair[4] = true;
            // Compute pair (m4, m4+1) corresponding to 4th solution.
            mpz_mul(m4, x, x);  // m4 = x^2 = w2^2
            mpz_mul(m4, m4, m2);            // m4 = m2*w2^2
            mpz_mul_ui(m4, m4, 4);          // m4 = 4*m2*w2^2
            mpz_out_str(outputfile, 10, m4);
            fputs(" 4", outputfile);
            fputs("\n", outputfile);
            // Check 8th solution:
            // corresponds to smooth pair if w8 = 8*x^4-8*x^2+1 = 8*x^2*(x^2-1)+1 is smooth.
            mpz_sub_ui(w8, x2, 1);      // w8 = x^2 - 1
            mpz_mul(w8, w8, x2);        // w8 = x^2*(x^2 - 1)
            mpz_mul_ui(w8, w8, 8);      // w8 = 8*x^2*(x^2 - 1)
            mpz_add_ui(w8, w8, 1);      // w8 = 8*x^2*(x^2 - 1) + 1
            if (is_smooth(w8, primes, num_primes)) {
                mpz_t m8;
                mpz_init(m8);

                is_pair[8] = true;
                // Compute (m8, m8+1).
                mpz_mul(m8, w4, w4);    // m8 = w4^2
                mpz_mul(m8, m8, m4);    // m8 = m4*w4^2
                mpz_mul_ui(m8, m8, 4);  // m8 = 4*m4*w4^2
                mpz_out_str(outputfile, 10, m8);
                fputs(" 8", outputfile);
                fputs("\n", outputfile);
                mpz_clear(m8);
            }
            mpz_clear(m4);
            mpz_clear(w8);
        }   
        mpz_clear(w4);     
    }
    
    // Check 3rd solution: 
    // corresponds to smooth pair if w30 = 2*x-1 and w31 = 2*x + 1 are smooth.
    mpz_add(w30, x, x);             // w30 = 2*x
    mpz_add_ui(w31, w30, 1);        // w31 = 2*x + 1
    mpz_sub_ui(w30, w30, 1);        // w30 = 2*x - 1
    if (is_smooth(w30, primes, num_primes) && is_smooth(w31, primes, num_primes)) {
        mpz_t m3, w90, w91;
        mpz_init(m3);
        mpz_init(w90);
        mpz_init(w91);

        is_pair[3] = true;
        // Compute pair (m3, m3+1).
        mpz_mul(m3, w31, w31);
        mpz_mul(m3, m3, m);     // m3 = m*w31^2
        mpz_out_str(outputfile, 10, m3);
        fputs(" 3", outputfile);
        fputs("\n", outputfile);
        // Check 9th solution:
        // corresponds to smooth pair if w90 = 8*x^3-6*x-1 = 2*x*(4*x^2-3)-1
        // and w91 = w90 + 2 are smooth.
        mpz_mul_ui(w90, x2, 4);
        mpz_sub_ui(w90, w90, 3);        // w90 = 4*x^2 - 3
        mpz_mul(w90, w90, x);     
        mpz_add(w90, w90, w90);         // w90 = 2*x*(4*x^2 - 3)
        mpz_sub_ui(w90, w90, 1);        // w90 = 2*x*(4*x^2 - 3) - 1
        mpz_add_ui(w91, w90, 2);        // w91 = 2*x*(4*x^2 - 3) + 1
        if (is_smooth(w90, primes, num_primes) && is_smooth(w91, primes, num_primes)) {
            mpz_t m9;
            mpz_init(m9);

            is_pair[9] = true;
            // Compute pair (m9, m9+1).
            mpz_mul(m9, w91, w91);
            mpz_mul(m9, m9, m3);     // m9 = m3*w91^2
            mpz_out_str(outputfile, 10, m9);
            fputs(" 9", outputfile);
            fputs("\n", outputfile);
            mpz_clear(m9);
        }
        // If also the second solution corresponds to a pair, 
        // check 6th solution.
        if (is_pair[2]) {
            mpz_t w6;
            mpz_init(w6);

            // Check 6th solution: 
            // corresponds to smooth pair if w6 = 4*x^2 - 3 is smooth.
            mpz_mul_ui(w6, x2, 4);
            mpz_sub_ui(w6, w6, 3);        // w6 = 4*x^2 - 3
            if (is_smooth(w6, primes, num_primes)) {
                mpz_t m6;
                mpz_init(m6);
                
                is_pair[6] = true;
                // Compute pair (m6, m6+1).
                mpz_mul(m6, w30, w31);
                mpz_mul(m6, m6, m6);
                mpz_mul(m6, m6, m2);    // m6 = m2*w30^2*w31^2
                mpz_out_str(outputfile, 10, m6);
                fputs(" 6", outputfile);
                fputs("\n", outputfile);

                // If the 6th and 4th solutions correspond to a smooth pair,
                // Check 12th solution.
                if (is_pair[4]) {
                    mpz_t w12;
                    mpz_init(w12);
                    // 12th sol corresponds to smooth pair if w12 = 16*x^4-16*x^2+1 
                    // = 16*x^2*(x^2-1)+1 is smooth.
                    mpz_sub_ui(w12, x2, 1);     // w12 = x^2 - 1
                    mpz_mul(w12, w12, x2);      // w12 = x^2*(x^2 - 1)
                    mpz_mul_ui(w12, w12, 16);   // w12 = 16*x^2*(x^2 - 1)
                    mpz_add_ui(w12, w12, 1);    // w12 = 16*x^2*(x^2 - 1) + 1
                    if (is_smooth(w12, primes, num_primes)) {
                        mpz_t m12;
                        mpz_init(m12);

                        is_pair[12] = true;
                        // Compute (m12, m12+1).
                        mpz_mul(m12, x, w6);
                        mpz_mul(m12, m12, m12);
                        mpz_mul(m12, m12, m6);
                        mpz_mul_ui(m12, m12, 4);    // m12 = 4*m6*x^2*w6^2
                        mpz_out_str(outputfile, 10, m12);
                        fputs(" 12", outputfile);
                        fputs("\n", outputfile);
                        mpz_clear(m12);
                    }
                    mpz_clear(w12);
                }
                mpz_clear(m6);
            }
            mpz_clear(w6);
        }
        mpz_clear(m3);
        mpz_clear(w90);
        mpz_clear(w91);
    }

    // Check 5th solution:
    // corresponds to smooth pair if w50 = 4*x^2-2*x-1 and w51 = 4*x^2+2*x-1
    // are smooth
    mpz_mul_ui(w51, x2, 4);     // w51 = 4*x^2
    mpz_sub_ui(w51, w51, 1);    // w51 = 4*x^2 - 1
    mpz_mul_ui(w50, x, 2);      // w50 = 2*x
    mpz_add(w51, w51, w50);     // w51 = 4*x^2 + 2*x - 1
    mpz_add(w50, w50, w50);     // w50 = 4*x
    mpz_sub(w50, w51, w50);     // w50 = 4*x^2 - 2*x - 1
    if (is_smooth(w50, primes, num_primes) && is_smooth(w51, primes, num_primes)) {
        mpz_t m5;
        mpz_init(m5);

        is_pair[5] = true;
        // Compute pair (m5, m5+1).
        mpz_mul(m5, w51, w51);
        mpz_mul(m5, m5, m);         // m5 = m*w51^2
        mpz_out_str(outputfile, 10, m5);
        fputs(" 5", outputfile);
        fputs("\n", outputfile);
        mpz_clear(m5);
        // If the 5th and 2nd solution correspond to smooth pairs
        // check the 10th solution.
        if (is_pair[2]) {
            mpz_t w10;
            mpz_init(w10);

            // Check 10th solution.
            // corresponds to a smooth pair if w10 = 16*x^4 - 20x^2 + 5 
            // = 4*x^2*(4*x^2 -5) + 5 is smooth.
            mpz_mul_ui(w10, x2, 4);     // w10 = 4*x^2
            mpz_sub_ui(w10, w10, 5);    // w10 = 4*x^2 - 5
            mpz_mul(w10, w10, x2);      // w10 = x^2*(4*x^2 - 5)
            mpz_mul_ui(w10, w10, 4);    // w10 = 4*x^2*(4*x^2 - 5)
            mpz_add_ui(w10, w10, 5);    // w10 = 4*x^2*(4*x^2 - 5) + 5
            if (is_smooth(w10, primes, num_primes)) {
                mpz_t m10;
                mpz_init(m10);

                is_pair[10] = true;
                // Compute (m10, m10+1).
                mpz_mul(m10, w50, w51);
                mpz_mul(m10, m10, m10);
                mpz_mul(m10, m10, m2);      // m10 = m2*w50^2*w51^2
                mpz_out_str(outputfile, 10, m10);
                fputs(" 10", outputfile);
                fputs("\n", outputfile);
                mpz_clear(m10);
            }
            mpz_clear(w10);
        }
    }

    // Check 7th solution:
    // corresponds to smooth pair if w70 = 8*x^3-4*x^2-4*x+1 
    // and w71 = 8*x^3+4*x^2-4*x-1 are smooth
    mpz_mul_ui(w70, x2, 4);
    mpz_sub_ui(w70, w70, 1);    // w70 = 4*x^2 - 1
    mpz_mul_ui(w71, x2, 2);
    mpz_sub_ui(w71, w71, 1);    // w71 = 2*x^2 -1
    mpz_mul(w71, w71, x);       // w71 = x*(2*x^2 -1)
    mpz_mul_ui(w71, w71, 4);    // w71 = 4*x*(2*x^2 -1) = 8*x^3 - 4*x
    mpz_add(w71, w71, w70);     // w71 = 8*x^3 + 4*x^2 - 4*x - 1
    mpz_add(w70, w70, w70);     // w70 = 2*(4*x^2 - 1)
    mpz_sub(w70, w71, w70);     // w70 = 8*x^3 - 4*x^2 - 4*x + 1
    if (is_smooth(w70, primes, num_primes) && is_smooth(w71, primes, num_primes)) {
        mpz_t m7;
        mpz_init(m7);

        is_pair[7];
        // Compute (m7, m7+1).
        mpz_mul(m7, w71, w71);
        mpz_mul(m7, m7, m);         // m7 = m*w71^2
        mpz_out_str(outputfile, 10, m7);
        fputs(" 7", outputfile);
        fputs("\n", outputfile);
        mpz_clear(m7);
    }

    // Check 11th solution:
    // corresponds to smooth pair if w110 = 32*x^5 - 16*x^4 - 32*x^3 + 12*x^2 + 6*x - 1
    // and w111 = 32*x^5 + 16*x^4 - 32*x^3 - 12*x^2 + 6*x + 1 are smooth.
    mpz_sub_ui(w110, x2, 1);        // w110 = x^2 - 1
    mpz_mul(w110, w110, x2);        // w110 = x^2*(x^2 - 1)
    mpz_mul_ui(w110, w110, 32);     // w110 = 32*x^2*(x^2 - 1)
    mpz_add_ui(w110, w110, 6);      // w110 = 32*x^2*(x^2 - 1) + 6
    mpz_mul(w110, w110, x);         // w110 = 32*x^3*(x^2 - 1) + 6*x = 32*x^5 - 32*x^3 + 6*x
    mpz_mul_ui(w111, x2, 4);       // w111 = 4*x^2
    mpz_sub_ui(w111, w111, 3);      // w111 = 4*x^2 - 3
    mpz_mul(w111, w111, x2);        // w111 = x^2*(4*x^2 - 3)
    mpz_mul_ui(w111, w111, 4);      // w111 = 4*x^2*(4*x^2 - 3)
    mpz_add_ui(w111, w111, 1);      // w111 = 4*x^2*(4*x^2 - 3) + 1 = 16*x^4 - 12*x^2 + 1
    mpz_sub(w110, w110, w111);      // w110 = 32*x^5 - 16*x^4 - 32*x^3 + 12*x^2 + 6*x - 1
    mpz_add(w111, w111, w111);      // w111 = 2*(16*x^4 - 12*x^2 + 1)
    mpz_add(w111, w110, w111);      // w111 = 32*x^5 + 16*x^4 - 32*x^3 - 12*x^2 + 6*x + 1
    if (is_smooth(w110, primes, num_primes) && is_smooth(w111, primes, num_primes)) {
        mpz_t m11;
        mpz_init(m11);

        is_pair[11] = true;
        // Compute (m11, m11+1).
        mpz_mul(m11, w111, w111);
        mpz_mul(m11, m11, m);       // m11 = m*w111^2
        mpz_out_str(outputfile, 10, m11);
        fputs(" 11", outputfile);
        fputs("\n", outputfile);
        mpz_clear(m11);
    }

    mpz_clear(w30);
    mpz_clear(w31);
    mpz_clear(w50);
    mpz_clear(w51);
    mpz_clear(w70);
    mpz_clear(w71);
    mpz_clear(w110);
    mpz_clear(w111);
    mpz_clear(m2);
    mpz_clear(x);   
    mpz_clear(x2);
}

bool check_higher_solutions(mpz_t m, mpz_t primes[], int num_primes, FILE *outputfile) {
    bool found = false;
    // Array of booleans for keeping track of whether 
    // consecutive solutions correspond to smooth pairs.
    bool is_pair[13];
    // The first two are not necessary, but included to have:
    // is_pair[i] = true if i-th solutions corresponds
    // to a smooth pair, else is_pair[i] = false.
    is_pair[0] = true;   
    is_pair[1] = true;
    for (int i=2; i<11; i++) {
        is_pair[i] = false;
    }
    // Initialize space for m2 as it might be used for
    // computing others such as m6, m10 and m12, not purely in the 
    // power-of-2 chain.
    mpz_t m2;
    mpz_init(m2);
    // Initialize space for wi0, wi1 for i in {3,5,7,11}. They are always 
    // computed to check for smoothness of pairs corresponding to the i-th
    // solutions.
    mpz_t w30, w31, w50, w51, w70, w71, w110, w111;
    mpz_init(w30);
    mpz_init(w31);
    mpz_init(w50);
    mpz_init(w51);
    mpz_init(w70);
    mpz_init(w71);
    mpz_init(w110);
    mpz_init(w111);
    
    // Compute x = 2*m+1 and its square x^2.
    mpz_t x, x2;
    mpz_init(x);
    mpz_init(x2);
    mpz_set_si(x, 1);
    mpz_addmul_ui(x, m, 2);
    mpz_mul(x2, x, x);
    
    // Now check smoothness of higher solutions by computing the wi.
    
    // Check 2nd solution:
    // corresponds to smooth pair if w2 = x is smooth.
    if (is_smooth(x, primes, num_primes)){
        mpz_t w4;
        mpz_init(w4);

        is_pair[2] = true;
        fputs(" 2", outputfile);
        found = true;
        
        // Check 4th solution:
        // corresponds to smooth pair if w4 = 2*x^2-1 is smooth.
        mpz_add(w4, x2, x2);
        mpz_sub_ui(w4, w4, 1);
        if (is_smooth(w4, primes, num_primes)){
            mpz_t w8;
            mpz_init(w8);

            is_pair[4] = true;
            fputs(" 4", outputfile);
            // Check 8th solution:
            // corresponds to smooth pair if w8 = 8*x^4-8*x^2+1 = 8*x^2*(x^2-1)+1 is smooth.
            mpz_sub_ui(w8, x2, 1);      // w8 = x^2 - 1
            mpz_mul(w8, w8, x2);        // w8 = x^2*(x^2 - 1)
            mpz_mul_ui(w8, w8, 8);      // w8 = 8*x^2*(x^2 - 1)
            mpz_add_ui(w8, w8, 1);      // w8 = 8*x^2*(x^2 - 1) + 1
            if (is_smooth(w8, primes, num_primes)) {
                is_pair[8] = true;
                fputs(" 8", outputfile);
            }
            mpz_clear(w8);
        }   
        mpz_clear(w4);     
    }
    
    // Check 3rd solution: 
    // corresponds to smooth pair if w30 = 2*x-1 and w31 = 2*x + 1 are smooth.
    mpz_add(w30, x, x);             // w30 = 2*x
    mpz_add_ui(w31, w30, 1);        // w31 = 2*x + 1
    mpz_sub_ui(w30, w30, 1);        // w30 = 2*x - 1
    if (is_smooth(w30, primes, num_primes) && is_smooth(w31, primes, num_primes)) {
        mpz_t w90, w91;
        mpz_init(w90);
        mpz_init(w91);

        is_pair[3] = true;
        fputs(" 3", outputfile);
        found = true;
        // Check 9th solution:
        // corresponds to smooth pair if w90 = 8*x^3-6*x-1 = 2*x*(4*x^2-3)-1
        // and w91 = w90 + 2 are smooth.
        mpz_mul_ui(w90, x2, 4);
        mpz_sub_ui(w90, w90, 3);        // w90 = 4*x^2 - 3
        mpz_mul(w90, w90, x);     
        mpz_add(w90, w90, w90);         // w90 = 2*x*(4*x^2 - 3)
        mpz_sub_ui(w90, w90, 1);        // w90 = 2*x*(4*x^2 - 3) - 1
        mpz_add_ui(w91, w90, 2);        // w91 = 2*x*(4*x^2 - 3) + 1
        if (is_smooth(w90, primes, num_primes) && is_smooth(w91, primes, num_primes)) {
            is_pair[9] = true;
            fputs(" 9", outputfile);
        }
        // If also the second solution corresponds to a pair, 
        // check 6th solution.
        if (is_pair[2]) {
            mpz_t w6;
            mpz_init(w6);

            // Check 6th solution: 
            // corresponds to smooth pair if w6 = 4*x^2 - 3 is smooth.
            mpz_mul_ui(w6, x2, 4);
            mpz_sub_ui(w6, w6, 3);        // w6 = 4*x^2 - 3
            if (is_smooth(w6, primes, num_primes)) {
                is_pair[6] = true;
                fputs(" 6", outputfile);

                // If the 6th and 4th solutions correspond to a smooth pair,
                // Check 12th solution.
                if (is_pair[4]) {
                    mpz_t w12;
                    mpz_init(w12);
                    // 12th sol corresponds to smooth pair if w12 = 16*x^4-16*x^2+1 
                    // = 16*x^2*(x^2-1)+1 is smooth.
                    mpz_sub_ui(w12, x2, 1);     // w12 = x^2 - 1
                    mpz_mul(w12, w12, x2);      // w12 = x^2*(x^2 - 1)
                    mpz_mul_ui(w12, w12, 16);   // w12 = 16*x^2*(x^2 - 1)
                    mpz_add_ui(w12, w12, 1);    // w12 = 16*x^2*(x^2 - 1) + 1
                    if (is_smooth(w12, primes, num_primes)) {
                        is_pair[12] = true;
                        fputs(" 12", outputfile);
                    }
                    mpz_clear(w12);
                }
            }
            mpz_clear(w6);
        }
        mpz_clear(w90);
        mpz_clear(w91);
    }

    // Check 5th solution:
    // corresponds to smooth pair if w50 = 4*x^2-2*x-1 and w51 = 4*x^2+2*x-1
    // are smooth
    mpz_mul_ui(w51, x2, 4);     // w51 = 4*x^2
    mpz_sub_ui(w51, w51, 1);    // w51 = 4*x^2 - 1
    mpz_mul_ui(w50, x, 2);      // w50 = 2*x
    mpz_add(w51, w51, w50);     // w51 = 4*x^2 + 2*x - 1
    mpz_add(w50, w50, w50);     // w50 = 4*x
    mpz_sub(w50, w51, w50);     // w50 = 4*x^2 - 2*x - 1
    if (is_smooth(w50, primes, num_primes) && is_smooth(w51, primes, num_primes)) {
        is_pair[5] = true;
        fputs(" 5", outputfile);
        found = true;
        // If the 5th and 2nd solution correspond to smooth pairs
        // check the 10th solution.
        if (is_pair[2]) {
            mpz_t w10;
            mpz_init(w10);

            // Check 10th solution.
            // corresponds to a smooth pair if w10 = 16*x^4 - 20x^2 + 5 
            // = 4*x^2*(4*x^2 -5) + 5 is smooth.
            mpz_mul_ui(w10, x2, 4);     // w10 = 4*x^2
            mpz_sub_ui(w10, w10, 5);    // w10 = 4*x^2 - 5
            mpz_mul(w10, w10, x2);      // w10 = x^2*(4*x^2 - 5)
            mpz_mul_ui(w10, w10, 4);    // w10 = 4*x^2*(4*x^2 - 5)
            mpz_add_ui(w10, w10, 5);    // w10 = 4*x^2*(4*x^2 - 5) + 5
            if (is_smooth(w10, primes, num_primes)) {
                is_pair[10] = true;
                // Compute (m10, m10+1).
                fputs(" 10", outputfile);
            }
            mpz_clear(w10);
        }
    }

    // Check 7th solution:
    // corresponds to smooth pair if w70 = 8*x^3-4*x^2-4*x+1 
    // and w71 = 8*x^3+4*x^2-4*x-1 are smooth
    mpz_mul_ui(w70, x2, 4);
    mpz_sub_ui(w70, w70, 1);    // w70 = 4*x^2 - 1
    mpz_mul_ui(w71, x2, 2);
    mpz_sub_ui(w71, w71, 1);    // w71 = 2*x^2 -1
    mpz_mul(w71, w71, x);       // w71 = x*(2*x^2 -1)
    mpz_mul_ui(w71, w71, 4);    // w71 = 4*x*(2*x^2 -1) = 8*x^3 - 4*x
    mpz_add(w71, w71, w70);     // w71 = 8*x^3 + 4*x^2 - 4*x - 1
    mpz_add(w70, w70, w70);     // w70 = 2*(4*x^2 - 1)
    mpz_sub(w70, w71, w70);     // w70 = 8*x^3 - 4*x^2 - 4*x + 1
    if (is_smooth(w70, primes, num_primes) && is_smooth(w71, primes, num_primes)) {
        is_pair[7];
        fputs(" 7", outputfile);
        found = true;
    }

    // Check 11th solution:
    // corresponds to smooth pair if w110 = 32*x^5 - 16*x^4 - 32*x^3 + 12*x^2 + 6*x - 1
    // and w111 = 32*x^5 + 16*x^4 - 32*x^3 - 12*x^2 + 6*x + 1 are smooth.
    mpz_sub_ui(w110, x2, 1);        // w110 = x^2 - 1
    mpz_mul(w110, w110, x2);        // w110 = x^2*(x^2 - 1)
    mpz_mul_ui(w110, w110, 32);     // w110 = 32*x^2*(x^2 - 1)
    mpz_add_ui(w110, w110, 6);      // w110 = 32*x^2*(x^2 - 1) + 6
    mpz_mul(w110, w110, x);         // w110 = 32*x^3*(x^2 - 1) + 6*x = 32*x^5 - 32*x^3 + 6*x
    mpz_mul_ui(w111, x2, 4);       // w111 = 4*x^2
    mpz_sub_ui(w111, w111, 3);      // w111 = 4*x^2 - 3
    mpz_mul(w111, w111, x2);        // w111 = x^2*(4*x^2 - 3)
    mpz_mul_ui(w111, w111, 4);      // w111 = 4*x^2*(4*x^2 - 3)
    mpz_add_ui(w111, w111, 1);      // w111 = 4*x^2*(4*x^2 - 3) + 1 = 16*x^4 - 12*x^2 + 1
    mpz_sub(w110, w110, w111);      // w110 = 32*x^5 - 16*x^4 - 32*x^3 + 12*x^2 + 6*x - 1
    mpz_add(w111, w111, w111);      // w111 = 2*(16*x^4 - 12*x^2 + 1)
    mpz_add(w111, w110, w111);      // w111 = 32*x^5 + 16*x^4 - 32*x^3 - 12*x^2 + 6*x + 1
    if (is_smooth(w110, primes, num_primes) && is_smooth(w111, primes, num_primes)) {
        is_pair[11] = true;
        // Compute (m11, m11+1).
        fputs(" 11", outputfile);
        found = true;
    }

    fputs("\n", outputfile);

    mpz_clear(w30);
    mpz_clear(w31);
    mpz_clear(w50);
    mpz_clear(w51);
    mpz_clear(w70);
    mpz_clear(w71);
    mpz_clear(w110);
    mpz_clear(w111);
    mpz_clear(m2);
    mpz_clear(x);   
    mpz_clear(x2);

    return found;
}