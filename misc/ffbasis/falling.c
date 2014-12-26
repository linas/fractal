/**
 * falling.c
 *
 * Falling factorial basis for the divisor function
 *
 * Linas Vepstas December 2014
 */

#include <math.h>
#include <stdlib.h>

#include "binomial.h"
#include "cache.h"
#include "falling.h"

/**
 * The operator that transforms binomial (falling factorial) into a
 * Dirichlet series.  Defined as:
 *    E_mk = (1/m - 1)^k
 *
 * We assume m>0 and that k >=0
 */
long double E_mk(unsigned int m, unsigned int k)
{
	long double em = m;
	em = (1.0L - em) / em;
	em = powl(em, ((long double) k));
	return em;
}

/**
 * A vector in the kernel of E aka a_k
 * If we did this correctly, a_k is given by solving
 *    sum_k=0^infty a_k x^k = sin (2pi/(1+x))
 * This returns a_k.
 * This breaks down around k=800 or so, not enough precision.
 */
long double a_k(unsigned int k)
{
	unsigned int n;
	long double sum = 0.0L;

	if (0 == k) return 0.0L;

	if (800 < k) fprintf(stderr, "Error: k=%d too large\n", k);

	DECLARE_LD_CACHE(a_kc);
	if (ld_one_d_cache_check(&a_kc, k))
		return ld_one_d_cache_fetch(&a_kc, k);

// WTF ?? -std=gnu should be enough??
#define M_PIl    3.141592653589793238462643383279502884L

	long double fourpi = - 4.0L * M_PIl * M_PIl;
	long double numer = - 2.0L * M_PIl;
	long double fact = 1.0L;
	for (n=0; n<3530; n++)
	{
		long double term = numer;
		long double bino = binomial (2*n+k, 2*n);
		bino /= fact;

		term *= bino;
		sum += term;
		// printf("n=%d numer=%g  bino=%g term=%g sum=%g\n", n, numer, bino, term, sum);

		if (10 < n && fabs(term/sum) < 1.0e-35) break;

		numer *= fourpi;
		fact *= (2.0L * n + 3.0L)  * (2.0L*n + 2.0L);
	}
	// printf(" ---------------------- above was k=%d\n", k);

	if (k%2 == 0) sum = - sum;

	ld_one_d_cache_store(&a_kc, sum, k);
	return sum;
}

