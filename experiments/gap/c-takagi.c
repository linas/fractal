/* 
 * c-takagi.C
 *
 * draw the Takagi curve for complex-valued w
 *
 * wow! checkout tagaki bumps!  wonder if one could make a sin(pi/x) 
 * approximation out of this...
 *
 * Linas October 2004
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

double bumps (double x)
{
	x -= floor (x);
	if (1.0 > x*3.0) return 0.0;
	if (2.0 < x*3.0) return 0.0;
	return triangle (3.0*x);
}

double saw (double x)
{
	x *= 3.0;
	if (1.0 > x) return 0.0;
	if (2.0 < x) return 0.0;
	return 2.0-x;
}

// #define Complex complex<long double>
#define Complex long double complex

/* The main, core basic takagi curve */
Complex takagi (Complex w, long double x)
{
	int k;
	Complex acc = 0.0L;
	Complex tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<36; k++)
	{
		Complex term = tw* triangle (tp*x);
		// long double term = tw* parabola_down (tp*x);
		// long double term = tw* parabola_up (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
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
	int nmax = 4723;
	// int nmax = 2048;

	if (argc <3)
	{
		printf ("Usage: %s <real w-value> <imag w-value>\n", argv[0]);
		exit (1);
	}
	double re_w = atof(argv[1]);
	double im_w = atof(argv[2]);

	Complex w = re_w + I* im_w;

	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);
		Complex tw = takagi (w, x);

#ifdef HALF_SYM
		// g
		Complex tw = takagi (w, 0.5*x);
		Complex ts = x + w*takagi (w, x);

		// g*g
		Complex tw = takagi (w, 0.25*x);
		Complex ts = 0.5*x + w*x + w*w*takagi (w, x);

		// g*g*g
		Complex tw = takagi (w, 0.125*x);
		Complex ts =  x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// r*g*g*g
		Complex tw = takagi (w, 1.0-0.125*x);
		Complex ts = x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// g^4
		Complex tw = takagi (w, 0.0625*x);
		Complex ts = 0.125*x*(1.0-16.0*w*w*w*w)/(1.0-2.0*w) + w*w*w*w*takagi (w, x);
		// gr
		Complex tw = takagi (w, 0.5-0.5*x);
		Complex ts = 1.0 - x + w*takagi (w,x);

		// grg
		Complex tw = takagi (w, 0.5-0.25*x);
		Complex ts = 1.0 +x*(-0.5+w) +w*w*takagi (w, x);

		// grg^3
		Complex tw = takagi (w, 0.5-0.0625*x);
		Complex ts = 1.0 + x*(-0.125+0.25*w +0.5*w*w + w*w*w) +w*w*w*w*takagi (w,x);

		double re_ts = creall (ts);
		double im_ts = cimagl (ts);
		
#endif
		double re_tw = creall (tw);
		double im_tw = cimagl (tw);
		
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, re_tw, im_tw);
		fflush (stdout);
	}
}
