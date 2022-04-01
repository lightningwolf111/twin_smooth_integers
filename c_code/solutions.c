#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "helpers.h"

void check_higher_solutions(mpz_t m, FILE *outputfile) {
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
    if (is_smooth(x)){
        mpz_t w4;
        mpz_init(w4);

        is_pair[2] = true;
        // Compute the pair (m2, m2+1) corresponding to the 2nd solution.
        mpz_add_ui(m2, m, 1);   // m2 = m+1
        mpz_mul(m2, m2, m);     // m2 = m*(m+1)
        mpz_mul_ui(m2, m2, 4);  // m2 = 4*m*(m+1)
        mpz_out_str(outputfile, 10, m2);
        // Check 4th solution:
        // corresponds to smooth pair if w4 = 2*x^2-1 is smooth.
        mpz_add(w4, x2, x2);
        mpz_sub_ui(w4, w4, 1);
        if (is_smooth(w4)){
            mpz_t m4, w8;
            mpz_init(w8);
            mpz_init(m4);

            is_pair[4] = true;
            // Compute pair (m4, m4+1) corresponding to 4th solution.
            mpz_mul(m4, x, x);  // m4 = x^2 = w2^2
            mpz_mul(m4, m4, m2);            // m4 = m2*w2^2
            mpz_mul_ui(m4, m4, 4);          // m4 = 4*m2*w2^2
            mpz_out_str(outputfile, 10, m4);
            // Check 8th solution:
            // corresponds to smooth pair if w8 = 8*x^4-8*x^2+1 = 8*x^2*(x^2-1)+1 is smooth.
            mpz_sub_ui(w8, x2, 1);      // w8 = x^2 - 1
            mpz_mul(w8, w8, x2);        // w8 = x^2*(x^2 - 1)
            mpz_mul_ui(w8, w8, 8);      // w8 = 8*x^2*(x^2 - 1)
            mpz_sub_ui(w8, w8, 1);      // w8 = 8*x^2*(x^2 - 1) + 1
            if (is_smooth(w8)) {
                mpz_t m8;
                mpz_init(m8);

                is_pair[8] = true;
                // Compute (m8, m8+1).
                mpz_mul(m8, w4, w4);    // m8 = w4^2
                mpz_mul(m8, m8, m4);    // m8 = m4*w4^2
                mpz_mul_ui(m8, m8, 4);  // m8 = 4*m4*w4^2
                mpz_out_str(outputfile, 10, m8);
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
    if (is_smooth(w30) && is_smooth(w31)) {
        mpz_t m3, w90, w91;
        mpz_init(m3);
        mpz_init(w90);
        mpz_init(w91);

        is_pair[3] = true;
        // Compute pair (m3, m3+1).
        mpz_mul(m3, w31, w31);
        mpz_mul(m3, m3, m);     // m3 = m*w31^2
        mpz_out_str(outputfile, 10, m3);
        // Check 9th solution:
        // corresponds to smooth pair if w90 = 8*x^3-6*x-1 = 2*x*(4*x^2-3)-1
        // and w91 = w90 + 2 are smooth.
        mpz_mul_ui(w90, x2, 4);
        mpz_sub_ui(w90, w90, 3);        // w90 = 4*x^2 - 3
        mpz_mul(w90, w90, x);     
        mpz_add(w90, w90, w90);         // w90 = 2*x*(4*x^2 - 3)
        mpz_sub_ui(w90, w90, 1);        // w90 = 2*x*(4*x^2 - 3) - 1
        mpz_add_ui(w91, w90, 2);        // w91 = 2*x*(4*x^2 - 3) + 1
        if (is_smooth(w90) && is_smooth(w91)) {
            mpz_t m9;
            mpz_init(m9);

            is_pair[9] = true;
            // Compute pair (m9, m9+1).
            mpz_mul(m9, w91, w91);
            mpz_mul(m9, m9, m3);     // m9 = m3*w91^2
            mpz_out_str(outputfile, 10, m9);
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
            if (is_smooth(w6)) {
                mpz_t m6;
                mpz_init(m6);
                
                is_pair[6] = true;
                // Compute pair (m6, m6+1).
                mpz_mul(m6, w30, w31);
                mpz_mul(m6, m6, m6);
                mpz_mul(m6, m6, m2);    // m6 = m2*w30^2*w31^2
                mpz_out_str(outputfile, 10, m6);

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
                    if (is_smooth(w12)) {
                        mpz_t m12;
                        mpz_init(m12);

                        is_pair[12] = true;
                        // Compute (m12, m12+1).
                        mpz_mul(m12, x, w6);
                        mpz_mul(m12, m12, m12);
                        mpz_mul(m12, m12, m6);
                        mpz_mul_ui(m12, m12, 4);    // m12 = 4*m6*x^2*w6^2
                        mpz_out_str(outputfile, 10, m12);
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
    if (is_smooth(w50) && is_smooth(w51)) {
        mpz_t m5;
        mpz_init(m5);

        is_pair[5] = true;
        // Compute pair (m5, m5+1).
        mpz_mul(m5, w51, w51);
        mpz_mul(m5, m5, m);         // m5 = m*w51^2
        mpz_out_str(outputfile, 10, m5);
        mpz_clear(m5);
        // If the 5th and 2nd solution correspond to smooth pairs
        // check the 10th solution.
        if (is_pair[2]) {
            mpz_t w10;
            mpz_init(w10);

            // Check 10th solution.
            // corresponds to a smooth pair if w10 = 16*x^4 - 20x^2 + 1 
            // = 4*x^2*(4*x^2 -5) + 1 is smooth.
            mpz_mul_ui(w10, x2, 4);     // w10 = 4*x^2
            mpz_sub_ui(w10, w10, 5);    // w10 = 4*x^2 - 5
            mpz_mul(w10, w10, x2);      // w10 = x^2*(4*x^2 - 5)
            mpz_mul_ui(w10, w10, 4);    // w10 = 4*x^2*(4*x^2 - 5)
            mpz_add_ui(w10, w10, 1);    // w10 = 4*x^2*(4*x^2 - 5) + 1
            if (is_smooth(w10)) {
                mpz_t m10;
                mpz_init(m10);

                is_pair[10] = true;
                // Compute (m10, m10+1).
                mpz_mul(m10, w50, w51);
                mpz_mul(m10, m10, m10);
                mpz_mul(m10, m10, m2);      // m10 = m2*w50^2*w51^2
                mpz_out_str(outputfile, 10, m10);
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
    if (is_smooth(w70) && is_smooth(w71)) {
        mpz_t m7;
        mpz_init(m7);

        is_pair[7];
        // Compute (m7, m7+1).
        mpz_mul(m7, w71, w71);
        mpz_mul(m7, m7, m);         // m7 = m*w71^2
        mpz_out_str(outputfile, 10, m7);
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
    mpz_mul_ui(w111, x2, 4)l;       // w111 = 4*x^2
    mpz_sub_ui(w111, w111, 3);      // w111 = 4*x^2 - 3
    mpz_mul(w111, w111, x2);        // w111 = x^2*(4*x^2 - 3)
    mpz_mul_ui(w111, w111, 4);      // w111 = 4*x^2*(4*x^2 - 3)
    mpz_add_ui(w111, w111, 1);      // w111 = 4*x^2*(4*x^2 - 3) + 1 = 16*x^4 - 12*x^2 + 1
    mpz_sub(w110, w110, w111);      // w110 = 32*x^5 - 16*x^4 - 32*x^3 + 12*x^2 + 6*x - 1
    mpz_add(w111, w111, w111);      // w111 = 2*(16*x^4 - 12*x^2 + 1)
    mpz_add(w111, w110, w111);      // w111 = 32*x^5 + 16*x^4 - 32*x^3 - 12*x^2 + 6*x + 1
    if (is_smooth(w110) && is_smooth(w111)) {
        mpz_t m11;
        mpz_init(m11);

        is_pair[11] = true;
        // Compute (m11, m11+1).
        mpz_mul(m11, w111, w111);
        mpz_mul(m11, m11, m);       // m11 = m*w111^2
        mpz_out_str(outputfile, 10, m11);
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