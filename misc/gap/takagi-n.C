/* 
 * takagi-n.C
 *
 * Draw non-dyadic variant of the Takagi curve
 *
 * Linas Decembner 2005
 */

#include <complex.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

long double triangle (long double x)
{
	long double t = x - floorl(x);
	if (0.5L > t) return 2.0L*t;
	return 2.0L*(1.0L-t);
}

long double balanced_triangle (long double x)
{
	return (triangle(x) - 0.5L);
}

double square (double x)
{
	x -= floor (x);
	if (1.0 > x*2.0) return -2.0;
	return 2.0;
}

long double parabola_down (long double x)
{
	long double t = x - floorl(x);
	return 4.0L*t*(1.0L-t);
}

long double parabola_up (long double x)
{
	long double t = x - floorl(x);
	if (0.5 < t) t = 1.0-t;
	return 4.0L*t*t;
}

/* Have triangle bumps at all integer frequencies */
long double n_takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	for (k=0; k<50; k++)
	{
		int d;
		long double term = tw* triangle (k+1,x);
		acc += term;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

// ===============================================
// ===============================================

main (int argc, char *argv[])
{
	int i;

	// int nmax = 512;
	// int nmax = 432;
	// int nmax = 1717;
	int nmax = 2048;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	double acc = 0.0;
	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);

		/* jitter, as this can make a difference */
		// x += ((double) rand()) / (RAND_MAX*((double)nmax));
		
		// double ts = isola (w, x);
		double tw = takagi (w, x);
		// double ts = tw;
		// tw  = exp (-tw);
		// acc += tw;
		// ts = acc;
		// double tw = iter_tak (w, x);

		// double tw = dirichlet_takagi (w, x);
		// double tw = plicative_takagi (w, x);
		// double ts = 2.0*gsl_sf_zeta (w) / 3.0 -0.5 - pow(2.0, -w)/3.0;

		double ts = sin_takagi (w, x);
		// double tw = dtakagi (w, x);
		// double tw = log (takagi(w,x));
		// double tw = ttakagi(w,x);
		// double tw = div_takagi (w, x);
		// double tw = takagi_bumps (w, x);
		// double tw = takagi_prime (w, x);

		// double tw = lytic (w, x);

#ifdef HALF_SYM
		// g
		double tw = takagi (w, 0.5*x);
		double ts = x + w*takagi (w, x);

		// g*g
		double tw = takagi (w, 0.25*x);
		double ts = 0.5*x + w*x + w*w*takagi (w, x);

		// g*g*g
		double tw = takagi (w, 0.125*x);
		double ts =  x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// r*g*g*g
		double tw = takagi (w, 1.0-0.125*x);
		double ts = x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// g^4
		double tw = takagi (w, 0.0625*x);
		double ts = 0.125*x*(1.0-16.0*w*w*w*w)/(1.0-2.0*w) + w*w*w*w*takagi (w, x);
		// gr
		double tw = takagi (w, 0.5-0.5*x);
		double ts = 1.0 - x + w*takagi (w,x);

		// grg
		double tw = takagi (w, 0.5-0.25*x);
		double ts = 1.0 +x*(-0.5+w) +w*w*takagi (w, x);

		// grg^3
		double tw = takagi (w, 0.5-0.0625*x);
		double ts = 1.0 + x*(-0.125+0.25*w +0.5*w*w + w*w*w) +w*w*w*w*takagi (w,x);

		double tw = takagi (w, 0.5*(1.0+x));
		double ts = 1.0-x + w*takagi (w, x);
		double tw = takagi (w, 1.0+0.0625*(x-1.0));
		double ts = 0.125*(1.0-x)*(1.0-16.0*w*w*w*w)/(1.0-2.0*w) + w*w*w*w*takagi (w, x);
		double tw = takagi (w, 0.75-0.25*x);
		double ts = 0.5 -2.0*w*w + x*(0.5-2.0*w +2.0*w*w);
		ts /= 1.0-2.0*w;
		ts += w*w*takagi (w, x);

		double tw = takagi_l (w, 1, 0.5*x);
		double ts = takagi_l (w, 0, 3.0*x/2.0);
		double tw = takagi_l (w, 2, 0.5*x);
		double ts = takagi_l (w, 1, 5.0*x/6.0);
		double tw = takagi_l (w, 3, 0.5*x);
		double ts = takagi_l (w, 2, 7.0*x/10.0);
		int k=5;
		double tw = takagi_l (w, k, 0.5*x);
		double ts = takagi_l (w, k-1, (2*k+1)*x/(2*(2*k-1)));

		int k=4;
		double tw = takagi_l (w, k, 0.5*x);
		double ts = (2*k+1)*x + w*takagi_l (w, k, x);
#endif
		// double tw = takagi_mink (w, x);
		// double ts = 1.0;

		// double tw = takagi (w, x);
		// double ts = takagi_prime (0.5*w, x);

		// double tw = takagi (w, 0.5-0.0625*x);
		// double ts = 0.125 + 0.25*w + 0.5*w*w - x*(0.125+0.25*w* +0.5*w*w - w*w*w);
		// ts += w*w*w*w*takagi (w, x);

		// double tw = takagi (w, 0.125*x);
		// double ts =  x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

#ifdef FOUR_DIM_REP
		// g
		double tw = takagi (w, 0.5*x);
		double ts = 2.0*x - x*x + w*takagi (w, x);

		// g*g
		double tw = takagi (w, 0.25*x);
		double ts = x*(1.0+2.0*w) - x*x*(w+0.25) + w*w*takagi (w, x);

		// g*r
		double tw = takagi (w, 0.5*(1.0-x));
		double ts = 1.0 -x*x + w*takagi (w, x);

		// r * g
		double tw = takagi (w, (1.0-0.5*x));
		double ts = 2.0*x -x*x + w*takagi (w, x);

		// grg^2
		double tw = takagi (w, (0.5-0.125*x));
		double ts = 1.0 + x*(w+2.0*w*w) - x*x*(w*w+0.25*w+0.0625) + w*w*w*takagi (w, x);

#endif
		// double tw = triangle (128.0*(0.5*x));
		// tw += triangle (128.0*(0.5+0.5*x));
		// tw *= 0.5;
		// double ts = triangle (64.0*x);

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
		fflush (stdout);
	}
}
