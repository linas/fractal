/*
 * binomial.c
 *
 * Miscellaneous binomial-coefficient related library routines.
 *
 * Linas Vepstas <linas@linas.org> Dec 2003, Dec 2004
 */


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"

// ======================================================
// brute-force factorial function
inline long double factorial (int n)
{
	int k;
	long double fac = 1.0L;
	for (k=2; k<=n; k++)
	{
		fac *= (long double) k;
	}
	if (0>n) fac = 0.0L;
	return fac;
}

// ======================================================
// brute-force binomial coefficent
// must have m<=n, m>=0
inline long double binomial (int n, int m)
{
	if (0>m) return 0.0L;
	if (2*m < n) m = n-m;
	int l = n-m;
	if (0>l) return 0.0L;
	int k;
	long double bin = 1.0L;
	for (k=1; k<=l; k++)
	{
		bin *= (long double) (m+k);
	}
	bin /= factorial (l);
	return bin;
}

#if TEST
void print_binomial (void)
{
	int i;
	for (i=0; i<10; i++)
	{
		int j;
		for(j=0; j<=i; j++)
		{
			printf ("%Lf ", binomial(i,j));
		}
		printf ("\n");
	}
	printf ("\n");
}
#endif

// ======================================================
// brute-force ratio of factorials
// return n! / m!  where m<=n

inline long double frat (int n, int m)
{
	if (m >= n) return 1.0L;
	
	int k;
	long double acc = 1.0;
	for (k=1; k<=n-m; k++)
	{
		acc *= (long double) (m+k);
	}
	return acc;
}

#if TEST
int test_frat (void)
{
	int i;
	int passtest = 1;
	for (i=0; i<45; i++)
	{
		int j;
		for(j=0; j<=i; j++)
		{
			int a = frat (i,j)+0.5;
			int b = factorial(i)/factorial(j) +0.5;
			if (a != b) { printf ("ERROR: test failed i=%d j=%d\n", i,j); passtest = 0; }
		}
	}
	if (passtest) printf ("frat test passed\n");
	printf ("\n");
	return !passtest;
}
#endif 

// ======================================================
// real-valued binomial coefficent
// must have m>=0
// returns z*(z-1)*(z-2)...*(z-m+1) / m!
//
inline long double fbinomial (long double z, int m)
{
	if(0>m) return 0.0L;

	int k;
	long double bin = 1.0L;
	long double fac = 1.0L;
	for (k=1; k<=m; k++)
	{
		bin *= z;
		// printf ("bin term k=%d bin=%Lg z=%Lg\n", k, bin, z);
		z -= 1.0L;
		fac *= (long double) k;

		// avoid exponent overflows with periodic divisions
		if (1.0e300 < fac)
		{
			bin /= fac;
			fac = 1.0L;
		}
	}
	bin /= fac;
	return bin;
}

// ======================================================
// complex-valued binomial coefficent
// must have m>=0
// returns z*(z-1)*(z-2)...*(z-m+1) / m!
//
inline long double complex cbinomial (long double complex z, int m)
{
	if(0>m) return 0.0L;

	int k;
	long double complex bin = 1.0L;
	long double fac = 1.0L;
	for (k=1; k<=m; k++)
	{
		bin *= z;
		// printf ("bin term k=%d bin=(%Lg,%Lg) z=(%Lg,%Lg)\n", 
		//    k, creall(bin), cimagl(bin), creall(z), cimagl(z));
		z -= 1.0L;
		fac *= (long double) k;

		// avoid exponent overflows with periodic divisions
		if (1.0e300 < fac)
		{
			bin /= fac;
			fac = 1.0L;
		}
	}
	bin /= fac;
	return bin;
}

#if TEST
int test_cbinomial (void)
{
	int n,k;

	int passtest = 1;
	for (n=0; n<380; n++)
	{
		for (k=0; k<=n; k++)
		{
			long double ba = binomial (n,k);
			long double complex ca = cbinomial (n,k);
			long double rca = creall (ca);
			long double ica = cimag (ca);
			if (fabsl (ica) > 1.0e-16)
			{
				printf ("Error: complex part not zero! n=%d k=%d\n",n,k);
				passtest = 0;
			}
			// printf ("n=%d k=%d ba=%Lg ba-ca=%Lg\n", n,k,ba, ba-ca);
			long double delt = rca - ba;
			delt /= ba;
			delt = fabsl(delt);
			if (1.0e-16 < delt)
			{
				printf ("Error: binomials dont match! n=%d k=%d ba=%Lg delt=%Lg\n",
									 n,k,ba,delt);
				passtest = 0;
			}
		}
	}
	if (passtest)
	{
		printf ("complex binomial test passed \n");
	}
	return !passtest;
}
#endif 


// ======================================================
// test harness
// 
#if TEST
main ()
{
	int failed_tests = 0;
	int total_tests = 0;
	failed_tests += test_frat ();  total_tests ++;
	failed_tests += test_cbinomial ();  total_tests ++;

	if (failed_tests)
	{
		printf ("Error: %d of %d tests failed\n", failed_tests, total_tests);
	}
	else
	{
		printf ("All %d tests passed\n", total_tests);
	}
}
#endif 

// ======================= END OF FILE ===================
