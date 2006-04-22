/* 
 * harmonic.h
 *
 * harmonic function, Riemann zeta and Hurwitz zeta function prototypes 
 */
#ifndef HARMONIC_H__
#define HARMONIC_H__

#ifdef __cplusplus
extern "C" {
#endif

// Euler's constant Abramowitz & Stegun Table 1.1 
#define M_GAMMA 0.577215664901532860606512L

#include <complex.h>

// Return Reiman zeta -1
// Must have integer n>=2 
long double zetam1 (int n);

/* Compute the Riemann zeta for general complex argument.
 * Uses the Hasse expansion. More or less accurate;
 * has some accuracy trouble near the pole at s=1;
 */
void riemann_zeta (double res, double ims, double *rez, double *imz);

// Return the trigamma function per A&S chapter on gamma
// trigamma = psi-prime derivative of psi
long double trigamma (long double x);

/* Return the n'th harmonic number */
long double harmonic (int n, long double ess);

/* Return Hurwitz zeta equiv of harmonic */
long double harmonic_hurwitz (int n, long double x, long double ess);

/* Return zeta function minus n'th harmonic */
long double zeta_minus_harmonic (int n, long double ess);

#ifdef __cplusplus
};
#endif

#endif /* HARMONIC_H__ */
