
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

long double swap12 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	r2 += a1;
	r2 = 1.0L / r2;
	r2 += a2;
	
	r2 = 1.0L / r2;
	return r2;
}

void grand (long double x, long double sre, long double sim, 
					 long double *pre, long double *pim)
{
#if 0
	long double ox = 1.0L/x;
	long double sw = ox - floorl(ox);
#else
	long double sw = swap12 (x);
#endif

	long double lnx = logl (x);

	long double ire = expl (sre*lnx);
	long double iim = ire;
	ire *= cosl (sim*lnx);
	iim *= sinl (sim*lnx);
	
	ire *= sw;
	iim *= sw;

	*pre = ire;
	*pim = iim;
}

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

	/* ??? huih ? */
	sum_re += 1.0L;

	*pre = sum_re;
	*pim = sum_im;
}

main ()
{
	long double sre, sim;
	long double zre, zim;

	sre = 0.5L;
	sim = 15.0L;

	int i;
	for (i=0; i<200; i++) {
		gral (120000, sre, sim, &zre, &zim);
		sim += 0.1;
		long double vv = zre*zre+zim*zim;
		printf ("%Lg\t%Lg\t%13.9Lg\t%13.9Lg\t%13.9Lg\n", sre, sim, zre, zim, vv);
	}
	
}

