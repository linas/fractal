
/* 
 * bern.c
 *
 * Try to look up bernoulli-based matrix elts
 * this is a flop, doesn't work
 * Linas Dec 2003
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "zetafn.h"



// ======================================================
// return the inverse bernoulli polynomial coefficient.
// return V_n,k where
// x^n = Sum(k=0 to n) V_n,k *B_k(x)
// 
// must have k<=n, k>=0
// 
double inverse_bernoulli (int n, int k)
{
	static int nmax = 0;
	static int ncache = -1;
	static double *arr = NULL;

	int idx = (n*(n+1))/2 + k;
	// do we need more memory for the cache ?
	if ((idx+n) > ncache) 
	{
		ncache = idx+n+500;
		arr = realloc(arr, ncache*sizeof (double));
		arr[0] = 1.0;
	}
	if (n > nmax)
	{
		int m;
		for (m=nmax; m<=n; m++)
		{
			int mdex = (m*(m+1))/2;
			arr[mdex+m] = 1.0;
			
			int j;
			for (j=0; j<m; j++)
			{
				int jdex = mdex + j;
				double acc = 0.0;
				int l;
				for (l=j; l<m; l++)
				{
					int ldex = (l*(l+1))/2 + j;
					acc += binomial (m,l) * bernoulli(m-l) * arr[ldex];
				}

				arr[jdex] = -acc;
			}
		}
		nmax = n;
	}
	return arr [idx];
}

#if TEST
main ()
{
	int i;
	for (i=0; i<15; i++)
	{
		int j;
		for(j=0; j<=i; j++)
		{
			printf ("%f ", inverse_bernoulli(i,j));
		}
		printf ("\n");
	}
	printf ("\n");
}
#endif


// return Reimann zeta for integer n>2
inline long double zeta (int n) { return 1.0L + zetam1(n); }

// ======================================================
// return matrix element of G

double gee (int m, int p)
{
	double val = frat (p+m+1, p+1) * zeta (p+m+2);
	val /= factorial (m);
	if (m%2==1) val = -val;
	return val;
}

double pee (int k, int l)
{
	int m;

	double sign = 1.0;
	if (k%2==1) sign = -1.0;
	
	double acc = 0.0;
	for (m=l; m<123; m++)
	{
		acc += inverse_bernoulli (m,l) * frat (m+k+1, m) *zeta (m+k+2);
		sign = -sign;
		printf (" peep m=%d acc=%g inv=%g\n", m, acc,
							 (m+1)*inverse_bernoulli (m,l)
							 );
	}
	return acc;
}

// ======================================================

static inline long double partial (int n, int m)
{
	long double acc = 0.0;
	int k;
	for (k=0; k<=n; k++)
	{
		acc += binomial (n,k) * bernoulli(n-k) * frat (m+k+1, k+1) *zeta (m+k+2);
	}
	return acc;
}

double ache (int l, int n)
{
	int m;

	long double sign = 1.0;
	if (l%2 == 1) sign = -1.0;

	long double acc = 0.0;
	for (m=l; m<73; m++)
	{
		acc += sign *inverse_bernoulli (m,l) * partial (n,m) / factorial (m);
		sign = -sign;
		printf ("m=%d acc=%g\n", m, (double) acc);
	}

	return acc;
}

// ======================================================

main ()
{
	// ache (6,6);
	//
	pee (0,0);
	printf ("=======\n");
	pee (0,1);
	return 0;
}

// ======================= END OF FILE ===================
