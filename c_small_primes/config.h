// number of primes below 32768
#define NUM_PRIMES 3512 
// smoothness bound
#define BOUND 32768
// The maximum period of continued fractions that we allow. 1000 should be very safe.
#define MAX_PERIOD 1000
// The cutoff for continued fractions. Will only find pairs up to this cutoff-1 bits from
// fundamental solutions.
#define BIT_CUTOFF (258)
// number of threads
#define NUM_THREADS 4