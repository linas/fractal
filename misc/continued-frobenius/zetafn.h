
/* zetafn.h
 *
 * function prototypes 
 */
#ifndef ZETAFN_H__
#define ZETAFN_H__

// Euler's constant Abramowitz & Stegun Table 1.1 
#define M_GAMMA 0.577215664901532860606512L;

#include <complex.h>

// return Reiman zeta -1
// must have n>=2 
long double zetam1 (int n);

// brute-force factorial function
long double factorial (int n);

// brute-force binomial coefficent 
// return n! / m! (n-m)!
// must have m<=n, m>=0
long double binomial (int n, int m);

// real-valued binomial coefficent
// must have m>=0
// returns z*(z-1)*(z-2)...*(z-m+1) / m!
//
long double fbinomial (long double z, int m);

// complex-valued binomial coefficent
// must have m>=0
// returns z*(z-1)*(z-2)...*(z-m+1) / m!
//
long double complex cbinomial (long double complex z, int m);

// brute-force ratio of factorials
// return n! / m!  where m<=n
long double frat (int n, int m);

// return the nth' bernoulli number
// this is a fairly fast algorithm ...
// must have n>=0
long double bernoulli (int n);

// return the trigamma function per A&S chapter on gamma
// trigamma = psi-prime derivative of psi
long double trigamma (long double x);

/* Return the n'th harmnic number */
long double harmonic (int n, long double ess);

/* Return Hurwitz zeta equiv of harmonic */
long double harmonic_hurwitz (int n, long double x, long double ess);

/* Return zeta function minus n'th harmonic */
long double zeta_minus_harmonic (int n, long double ess);

#endif /* ZETAFN_H__ */
