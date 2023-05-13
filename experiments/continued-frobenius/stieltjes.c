
/*
 * Verify alternate derivation of the stieltjes constants
 *
 * Linas Sept 2004
 */

#include "ache.h"
#include "zetafn.h"
#include <math.h>
#include <stdio.h>

long double
stieltjes_gamma (int k)
{
	int n;
	long double sgn;
	long double gam = 0.0;
	if (0>k) return 0.0;
	if (0==k) return M_GAMMA;

	gam = a_sub_n(k-1)*((double) k);
	printf ("first term=%Lg\n", gam);

	sgn = 1.0;
	if (k%2==1) sgn = -1.0;
sgn=-1.0;
	for (n=0; n<30; n++)
	{
		long double term;
		term = stirling_first (n+k,k) + stirling_first(n+k,k-1);
		printf ("term n=%d stir=%Lg   ", n, term);
		term /= frat(n+k, k);
		printf ("w/fact=%Lg   ", term);
		term *= a_sub_n(n+k);
		printf ("term==%Lg\n", sgn*term);
		gam += sgn*term;
		sgn = -sgn;
	}
	return gam;
}

int
main ()
{
	double g;
	int k;

	k = 8;
	g = stieltjes_gamma (k);
	printf ("Stieltjes (%d)=%20.15g\n", k, g);
	
	return 0;
}

