/*
 * harmonic.c
 *
 * Miscellaneous zeta-function and related library routines.
 *
 * Linas Vepstas <linas@linas.org> Dec 2003, Dec 2004
 */


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

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

/* Return the n'th harmnic number */
long double harmonic (int n, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = 1; k<= n; k++)
	{
		long double kay = k;
		acc += powl (kay, -ess);
	}
	return acc;
}

/* Return Hurwitz zeta equiv of harmonic */
long double harmonic_hurwitz (int n, long double x, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = 1; k<= n; k++)
	{
		long double kay = k;
		kay += x;
		acc += powl (kay, -ess);
	}
	return acc;
}

/* Return zeta function minus n'th harmonic */
long double zeta_minus_harmonic (int n, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = n+1; k<10000; k++)
	{
		long double kay = k;
		long double term = powl (kay, -ess);
		acc += term;
		if (term < 1.0e-20*acc) break;
	}
	return acc;
}

// ======================================================
/* Compute the Riemann zeta for general complex argument.
 * Uses the Hasse expansion. More or less accurate;
 * has some accuracy trouble near the pole at s=1;
 */
void riemann_zeta (double res, double ims, double *rez, double *imz)
{
	double zre = 0.0;
	double zim = 0.0;
	double err;

	int n;

	double twn = 0.5;

	/* don't change the 110 -- it seems liek the right 
	 * thing for stuff going in the imaginary direction
	 */
	for (n=0; n<110; n++)
	{
		int k;
		double reb = 0.0;
		double imb = 0.0;
		double sgn = 1.0;
		for (k=0; k<=n; k++)
		{
			double r = binomial (n,k);
			// printf ("duude %d %d bin=%g\n", n, k, r);
			double lnk = log (k+1.0);
			r *= sgn * exp (-res*lnk);
			reb += r * cos (ims*lnk);
			imb -= r * sin (ims*lnk);
			sgn = -sgn;
		}

		double ret = twn * reb;
		double imt = twn * imb;
		zre += ret;
		zim += imt;

		err = ret*ret+imt*imt;
		// printf ("duude n=%d z=(%g %g) err=%g\n", n, zre, zim, err);

		/* Along imaginary axies, the error never seems 
		 * get less than 1.0e-34, no matter what (I don't get why)
		 * so cut off here.
		 */
		if (err < 1.0e-33) break;
		twn *= 0.5;
	}

	// printf ("duude leave with n=%d err=%g \n", n, sqrt (err));

	double r = pow (2.0, 1-res);
	double ret = 1.0 - r* cos (ims*M_LN2);
	double imt = r* sin (ims*M_LN2);
	r = ret*ret + imt*imt;
	r = 1.0/r;

	ret *= r;
	imt = -imt*r; 

	double tmp = ret * zre - imt *zim;
	zim = ret*zim + imt *zre;
	zre = tmp;
	
	*rez = zre;
	*imz = zim;
}

#ifdef TEST_ZETA

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

int
main (int argc, char * argv[])
{
	double s,t=0.0;
	int n;

	int error_occured = 0;

	for (n=2; n<=40; n++)
	{
		double reg, img;
		riemann_zeta (n, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta_int (n);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-15) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for n=%d   error=%g %g \n", n, err,  img);
			error_occured ++;
		}
	}

	for (t=2.0; t<=66.0; t+=0.0314683)
	{
		double reg, img;
		riemann_zeta (t, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta (t);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-15) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for s=%g   error=%g %g \n", t, err,  img);
			error_occured ++;
		}
	}

	for (t=1.01; t<=2.1; t+=0.00314683)
	{
		double reg, img;
		riemann_zeta (t, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta (t);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-13) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for s=%g   error=%g %g \n", t, err,  img);
			error_occured ++;
		}
	}

	for (s=-40.0; s<=40.0; s += 0.4356346)
	{
		for (t=0.0; t<=48.0; t+=0.6314683)
		{
			double reg, img;
			riemann_zeta (0.5, t, &reg, &img);
			double nreg, nimg;
			riemann_zeta (0.5, -t, &nreg, &nimg);
	
			double rerr = reg-nreg;
			double ierr = img+nimg;
			if ((fabs(rerr) > 1.0e-13) || (fabs (ierr) > 1.0e-15))
			{
				printf ("ERROR for s=%g   error=%g %g \n", t, rerr,  ierr);
				error_occured ++;
			}
		}
	}

	return error_occured;
}
#endif /* TEST_ZETA */

// ======================================================
// test harness
// 
#if TEST
main ()
{
	int failed_tests = 0;
	int total_tests = 0;
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

