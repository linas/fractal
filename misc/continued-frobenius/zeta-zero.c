
/*
 * zeta-sero.c
 *
 * explore zeros of generalized zeta function
 * generalized zeta function equals zeta when z=1/2
 * (currently borken)
 *
 * Linas Dec 2003
 */

#include <complex.h>
#include <math.h>

#include "zetafn.h"

#define PC(k,S,C)  \
	printf ("%d %s %s = %Lg +i %Lg)\n", k, S, #C, creall(C), cimagl(C));

// generalized zeta function equals zeta when z=1/2
long double complex geta (long double complex s, long double complex z)
{
	int k;

	long double complex gacc = 0.0L;

	long double complex zeke = z;
	long double sign = -1.0L;
	for (k=1; k<130; k++)
	{
		long double complex term = zetam1 (k+1);
		PC(k,"zeta", term);
		term *= (long double) k;
		term = 1.0L - term;
		term /= (long double) (k*(k+1));

		// xxxxxxxxxxxxxx
		term = zetam1 (k+1);
		term /= (long double) (k);
		
		PC(k,"zeta-kk", term);
		long double complex cbin = cbinomial (s-1.0L, k);
		PC(k,"cbin", cbin);
		term *= cbin;
		PC(k,"binom zeta-kk", term);
		term *= zeke;
		PC(k,"zeke binom zeta-kk", term);
		term *= sign;

		PC(k,"sign", term);
		gacc += term;
		PC(k,"acc", gacc);
		printf ("----------\n");
		
		zeke *= z;
		sign = -sign;
	}
	gacc += 1.0L - M_GAMMA;
	gacc *= - cpow (z, 1.0L-s);
	gacc += 1.0L / (s-1.0L);
	gacc *= s;
	return gacc;
}

main ()
{
	long double complex ess = 0.5L + 14.1347IL;
	ess = 2.01L;
	long double complex zee = 0.5L;
	PC (0, "yo", ess);
	
	long double complex g = geta (ess, zee);

	printf ("its %Lg %Lg\n", creall (g), cimagl(g));
}
