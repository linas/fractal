
/*
 * swap.c
 *
 * Integral the permuation group of continued fractions
 * Expect to get riemann zeta in the gauss map case
 * and that is what we seem to get ... need high integration 
 * order though to get anything on the r=1/2 axis ... 
 *
 * Linas Feb 2005
 */ 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

inline long double swap12 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	r2 += a1;
	r1 = 1.0L / r2;
	r1 += a2;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap13 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-15 > r2) return 0.0;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;

	r3 += a1;

	r2 = 1.0L / r3;
	r2 += a2;
	
	r1 = 1.0L / r2;
	r1 += a3;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap23 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-15 > r2) return 0.0;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;

	r3 += a2;

	r2 = 1.0L / r3;
	r2 += a3;
	
	r1 = 1.0L / r2;
	r1 += a1;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap14 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-15 > r2) return 0.0;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;
	if (1.0e-15 > r2) return 0.0;
	ox = 1.0L/r3;
	long double a4 = floorl(ox);
	long double r4 = ox - a4;

	r4 += a1;

	r3 = 1.0L / r4;
	r3 += a3;
	
	r2 = 1.0L / r3;
	r2 += a2;
	
	r1 = 1.0L / r2;
	r1 += a4;
	
	x = 1.0L / r1;
	return x;
}

/* The integrand, which is swap(x) * x^s */
void grand (long double x, long double sre, long double sim, 
					 long double *pre, long double *pim)
{

	long double lnx = logl (x);

	long double ire = expl (sre*lnx);
	long double iim = ire;
	long double phi = sim*lnx;
	ire *= cosl (phi);
	iim *= sinl (phi);
	
#if 0
	long double ox = 1.0L/x;
	long double sw = ox - floorl(ox);
#else
	long double sw = swap12 (x);
#endif
	ire *= sw;
	iim *= sw;

	*pre = ire;
	*pim = iim;
}

/* compute multiple integrands at once for performance */
void multi_grand (long double x, long double sre, long double sim, 
					 long double pre[5], long double pim[5])
{

	long double lnx = logl (x);

	long double ire = expl (sre*lnx);
	long double iim = ire;
	long double phi = sim*lnx;
	ire *= cosl (phi);
	iim *= sinl (phi);
	
	long double ox = 1.0L/x;
	long double sw = ox - floorl(ox);

	pre[0] = sw*ire;
	pim[0] = sw*iim;

	sw = swap12 (x);
	pre[1] = sw*ire;
	pim[1] = sw*iim;

	sw = swap13 (x);
	pre[2] = sw*ire;
	pim[2] = sw*iim;

	sw = swap23 (x);
	pre[3] = sw*ire;
	pim[3] = sw*iim;

	sw = swap14 (x);
	pre[4] = sw*ire;
	pim[4] = sw*iim;
}

/* compute single integral of the integrand.
 * actually compute 
 * zeta = s/(s-1) - s \int_0^1 swap(x) x^{s-1} dx 
 */
void gral(int nsteps, long double sre, long double sim, 
					 long double *pre, long double *pim)
{
	int i;

	long double step = 1.0L / ((long double) nsteps);
	
	long double sum_re= 0.0L;
	long double sum_im= 0.0L;
	long double x = 1.0L - 0.5*step;

	/* integrate in a simple fashion */
	for (i=0; i<nsteps; i++)
	{
		long double val_re, val_im;
		grand (x, sre-1.0, sim, &val_re, &val_im);
		sum_re += val_re;
		sum_im += val_im;
		x -= step;
	}
	sum_re *= step;
	sum_im *= step;

	/* mult by s */
	long double tmp = sum_re * sre - sum_im*sim;
	sum_im = sum_im *sre + sum_re *sim;
	sum_re = tmp;

	/* compute 1/(s-1) */
	sre -= 1.0L;
	long double v = sre*sre+sim*sim;
	sre /= v;
	sim = -sim/v;

	/* subtract */
	sum_re = sre - sum_re;
	sum_im = sim - sum_im;

	/* s/(s-1) so add 1 now */
	sum_re += 1.0L;

	*pre = sum_re;
	*pim = sum_im;
}

/* compute multiple integrals */
void multi_gral(int nsteps, long double sre, long double sim, 
					 long double sum_re[5], long double sum_im[5])
{
	int i, j;
	long double pre[5];
	long double pim[5];

	long double step = 1.0L / ((long double) nsteps);
	
	for (j=0; j<5; j++) 
	{
		sum_re[j]= 0.0L;
		sum_im[j]= 0.0L;
	}

	long double r = RAND_MAX;
	r = 1.0L / r;
	/* Integrate in a simple fashion */
	for (i=0; i<nsteps; i++)
	{
		int nr = rand();
		if (0 == nr) continue;
		long double x = (long double) nr;
		x *= r;

		multi_grand (x, sre-1.0, sim, pre, pim);
		for (j=0; j<5; j++) 
		{
			sum_re[j] += pre[j];
			sum_im[j] += pim[j];
		}
	}

	/* compute 1/(s-1) */
	double osre = sre - 1.0L;
	long double v = osre*osre+sim*sim;
	osre /= v;
	double osim = -sim/v;
	
	for (j=0; j<5; j++) 
	{
		sum_re[j] *= step;
		sum_im[j] *= step;

		/* mult by s */
		long double tmp = sum_re[j] * sre - sum_im[j]*sim;
		sum_im[j] = sum_im[j] *sre + sum_re[j] *sim;
		sum_re[j] = tmp;
	
		/* subtract */
		sum_re[j] = osre - sum_re[j];
		sum_im[j] = osim - sum_im[j];

		/* s/(s-1) so add 1 now */
		sum_re[j] += 1.0L;
	}
}

#ifdef SINGLE_MAIN
int
main (int argc, char * argv[])
{
	long double sre, sim;
	long double zre, zim;

	sre = 0.5L;
	sim = 13.0L;

	int i;
	for (i=0; i<200; i++) {
		gral (204123, sre, sim, &zre, &zim);
		sim += 0.1;
		long double vv = zre*zre+zim*zim;
		printf ("%Lg\t%Lg\t%13.9Lg\t%13.9Lg\t%13.9Lg\n", sre, sim, zre, zim, vv);
		fflush (stdout);
	}
	return 0;
}

#else 

int
main (int argc, char * argv[])
{
	long double sre, sim;
	long double zre[5], zim[5];

	sre = 0.5L;
	sim = 4.0L;

	int npts = atoi (argv[1]);

	printf ("#\n# col 1-2: re s,im s \n");
	printf ("#\n# col 3-4: re & im  h \n");
	printf ("#\n# col 5-6: re & im  swap 12\n");
	printf ("#\n# col 7-8: re & im  swap 13\n");
	printf ("#\n# col 9-10: re & im  swap 23\n");
	printf ("#\n# col 11-12: re & im  swap 14\n");
	printf ("#\n# RANDMAX = %d\n", RAND_MAX);
	printf ("#\n# n integration samples =%d\n", npts);
	printf ("#\n#\n");
	int i, j;
	for (i=0; i<300; i++) {
		multi_gral (npts, sre, sim, zre, zim);
		sim += 0.1;
		printf ("%Lg\t%Lg", sre, sim);
		for (j=0; j<5; j++) {
			printf ("\t%13.9Lg\t%13.9Lg", zre[j], zim[j]);
		}
		printf ("\n");
		fflush (stdout);
	}
	
	return 0;
}
#endif

