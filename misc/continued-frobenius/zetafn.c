
/*
 * zetafn.c
 *
 * Miscellaneous zeta-function and related library routines.
 *
 * Linas Vepstas <linas@linas.org> Dec 2003
 */


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "zetafn.h"

/* zeta function for integer values */
#define NZ 43
static long double zeta[NZ];

static void init_zeta (void)
{
	/* zeta minus 1 copied from abramowitz & stegun */
	zeta[2] = 0.64493406684822643647;
	zeta[3] = 0.20205690315959428540;
	zeta[4] = 0.08232323371113819152;
	zeta[5] = 0.03692775514336992633;
	zeta[6] = 0.01734306198444913971;
	
	zeta[7] = 0.00834927738192282684;
	zeta[8] = 0.00407735619794433938;
	zeta[9] = 0.00200839282608221442;
	zeta[10] = 0.00099457512781808534;
	zeta[11] = 0.00049418860411946456;
	zeta[12] = 0.00024608655330804830;
	
	zeta[13] = 0.00012271334757848915;
	zeta[14] = 0.00006124813505870483;
	zeta[15] = 0.00003058823630702049;
	zeta[16] = 0.00001528225940865187;
	zeta[17] = 0.00000763719763789976;
	zeta[18] = 0.00000381729326499984;
	
	zeta[19] = 0.00000190821271655394;
	zeta[20] = 0.00000095396203387280;
	zeta[21] = 0.00000047693298678781;
	zeta[22] = 0.00000023845050272773;
	zeta[23] = 0.00000011921992596531;
	zeta[24] = 0.00000005960818905126;
	
	zeta[25] = 0.00000002980350351465;
	zeta[26] = 0.00000001490155482837;
	zeta[27] = 0.00000000745071178984;
	zeta[28] = 0.00000000372533402479;
	zeta[29] = 0.00000000186265972351;
	zeta[30] = 0.00000000093132743242;

	zeta[31] = 0.00000000046566290650;
	zeta[32] = 0.00000000023283118337;
	zeta[33] = 0.00000000011641550173;
	zeta[34] = 0.00000000005820772088;
	zeta[35] = 0.00000000002910385044;
	zeta[36] = 0.00000000001455192819;

	zeta[37] = 0.00000000000727595984;
	zeta[38] = 0.00000000000363797955;
	zeta[39] = 0.00000000000181898965;
	zeta[40] = 0.00000000000090949478;
	zeta[41] = 0.00000000000045474738;
	zeta[42] = 0.00000000000022737368;
}

static int zeta_not_init = 1;

/* return Reimann zeta(n) -1 */
long double zetam1 (int n)
{
	if (zeta_not_init) 
	{
		zeta_not_init = 0;
		init_zeta();
	}
	if (n<2)
	{
		printf ("Error bad zeta %d\n",n);
		return 0.0L;
	}
	if (n<23)
	{
		return zeta[n];
	}

	long double twop = 1.0;
	long double threp = 1.0;
	long double fourp = 1.0;
	long double fivep = 1.0;
	long double sixp = 1.0;
	long double sevp = 1.0;
	long double eigp = 1.0;
	
	int i;
	for (i=0; i<n; i++)
	{
		twop *= 0.5;
		threp *= 1.0/3.0;
		fourp *= 0.25;
		fivep *= 0.2;
		sixp *= 1.0/6.0;
		sevp *= 1.0/7.0;
		eigp *= 0.125;
	}
	long double rv = twop + threp + fourp + fivep +sixp + sevp + eigp;
	// printf ("n=%d dif=%Lg per=%Lg\n", n, rv-zeta[n], (rv-zeta[n])/rv);
	return rv;
}

// =================================================================

#if TESTING
main ()
{
	int i;
	for (i=2; i<43; i++)
	{
		zetam1 (i);
	}
}
#endif

#ifdef GSL_COMPARE_TEST

#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_zeta.h>
#include "zetafn.h"

main ()
{
	int i;

	for (i=2; i<70; i++)
	{
		int rc;
		gsl_sf_result res;
		rc = gsl_sf_zeta_int_e (i, &res);

		long double zm1;
		zm1 = zetam1 (i);
		
		long double gm1 = res.val;
		gm1 -= 1.0L;

		long double diff = zm1 - gm1;
		printf ("i=%d zm1=%Lg diff=%Lg rat=%Lg\n", i, zm1, diff, diff/zm1);
	}
}
#endif

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
// return the nth' bernoulli number
// this is a fairly fast algorithm ...
// must have n>=0
long double bernoulli (int n)
{
	static int nmax = 0;
	static int ncache = 0;
	static long double *bern_cache = NULL; // cache of precomputed values

	if (1 > n) return 1.0;
	if (1 == n) return -0.5;
	if (n%2 == 1) return 0.0;
	
	if (n > ncache)  // do we need more memory for the cache?
	{
		ncache = 2*n+500;
		bern_cache = realloc(bern_cache, ncache*sizeof (long double));
	}

	// is the requested value higher than what we have in cache?
	if (n > nmax) 
	{
		if (0 == nmax)
		{
			bern_cache[0] = 1.0L;
			bern_cache[1] = -0.5L;
			bern_cache[2] = 1.0L/6.0L;
			nmax = 2;
		}
		int k;

		// use recursion relation on bernoulli numbers
		for (k=nmax+2; k<n+1; k+=2) 
		{
			long double acc = 0.0L;
			acc =  bern_cache[0];
			acc += ((long double) (k+1)) * bern_cache[1];
			acc += 0.5L * ((long double)(k+1)*k) * bern_cache[2];
			int j;
			for (j=4; j<k-1; j+=2)
			{
				acc += binomial (k+1, j) * bern_cache[j];
			}
			acc /= (long double) (k+1);
			bern_cache[k] = -acc;
		}

		nmax = n; // note nmax always even here
	}

	return bern_cache[n];
}

#if TEST
int test_bernoulli (void)
{
	int i;
#define L_PI 3.1415926535897932384626433832795029L
	long double lp = 2.0L * L_PI;
	lp = logl (lp);
	long double l2 = logl (2.0L);
	int passtest = 1;
	for (i=2; i<150; i+=2)
	{
		long double bound = 2*factorial(i) / expl (((long double) i) * lp);
		long double upper = expl (((long double)(1-i)) * l2);
		upper = 1.0L/(1.0L-upper);
		upper *= bound;
		long double bern = bernoulli (i);
		if (i<50) {
			printf ("%d bernoulli=%Lf upper bound = %Lf lower bound=%Lf \n", 
							 i, bern, upper, bound);
		}
		if (i%4 == 0) bern = -bern;
		if (upper <= bern) 
		{
			printf ("ERROR: fail upper bound at i=%d\n", i);
			passtest = 0;
		}
		if (bound >= bern) 
		{
			printf ("ERROR: fail lower bound at i=%d\n", i);
			passtest = 0;
		}
	}
	if (passtest) printf ("Bernoulli test passed \n");
	printf ("\n");

	return !passtest;
}
#endif 

// ======================================================
// trigamma function is the first derivative of digamma, per A&S

long double 
trigamma (long double x)
{
	long double y = 2.0L - x;
	long double fly = floorl (y);
	y -= fly;
	if (0.5 < y) { y -= 1.0L; fly += 1.0L; }

	// For this loop to converge, we want -2 <y < +2
	// The best convergence is for -0.5<y<+0.5 
	int n;
	long double acc = 0.0;
	long double ym = 1.0L;
	for (n=0; n<50; n++)
	{
		long double term = (n+1);
		term *= zetam1(n+2) * ym;
		acc += term;
		ym *= y;
	}

	// Now go back: psi(1+z) = -1/z^2 + psi (z)
	n = fly;
	if (0<n)
	{
		int k;
		long double ex = 2.0L - y;
		for (k=0; k<n; k++)
		{
			long double term = -(k+1);
			term += ex;
			acc += 1.0L / (term*term);
		}
	}
	else if (0 > n)
	{
		n = -n;
		int k;
		long double ex = 2.0L - y;
		for (k=0; k<n; k++)
		{
			long double term = k;
			term += ex;
			acc -= 1.0L / (term*term);
		}
	}

	return acc;
}

#if TEST

int test_trigamma (void)
{
	int passtest = 1;
	
	// trigamma should obey the multiplication formula, A&S eqn 6.4.8  
	// so test against that for a variety of inputs
	long double z;
	for (z=-0.95; z<3.0; z+=0.1)
	{
		for (m=2; m<15; m++)
		{
			long double acc = 0.0L;
			long double emm = 1.0/ (long double) m;
			int k;
			for (k=0; k<m; k++)
			{
				long double term = z + k*emm;
				term = trigamma (term);
				acc += term;
			}
			acc *= emm*emm;
			long double comp = trigamma (((long double) m) *z);
			acc -= comp;
			acc /= comp;
			if (2.0e-15 < fabs (acc)) 
			{
				printf ("Error: trigamma test fails at m=%d z=%Lg per=%Lg\n", m, z, acc);
				passtest = 0;
			}
		}
	}
	
	if (passtest) printf ("Trigamma test passed \n");
	printf ("\n");

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
	failed_tests += test_bernoulli ();  total_tests ++;
	failed_tests += test_trigamma ();  total_tests ++;

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
