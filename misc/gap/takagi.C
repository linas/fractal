
/* 
 * takagi.C
 *
 * draw the Takagi curve
 *
 * wow! checkout tagaki bumps!  wonder if one could make a sin(pi/x) 
 * approximation out of this...
 *
 * Linas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double triangle (double x)
{
	double t = x - floor(x);
	if (0.5 > t) return 2.0*t;
	return 2.0*(1.0-t);
}

double square (double x)
{
	x -= floor (x);
	if (1.0 > x*2.0) return 0.0;
	return 1.0;
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

double takagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000000; k++)
	{
		acc += tw* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

double sin_takagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<10000000; k++)
	{
		acc += tw* (1.0 - cos(tp*x*(2.0*M_PI)))*0.5;
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

double dtakagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 4.0;
	for (k=2; k<10000000; k++)
	{
		double term = k* (k-1)*tw;
		acc += term* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > term) break;
	}

	return acc;
}

double ttakagi (double s, double x)
{
	int k;
	double exps = exp(-s);
	double w = exps;
printf ("# duude w=%g\n", w);
	double acc = 0.0;
	double tw = 1.0;
	double tp = 2.0;
	for (k=1; k<10000000; k++)
	{
		acc += k* tw* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc*exps;
}

double div_takagi (double cut, double x)
{
	int k;
	double w = 2.0;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<100000000; k++)
	{
		double reg = exp (-k*k*cut*cut);
		double term = reg * tw;
		acc += term * triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > term) break;
	}

	acc *= 2.0 *cut*cut;

	return acc;
}

double takagi_bumps (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* bumps (tp*x);
		// acc += tw* saw (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

double takagi_prime (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* square (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

double takagi_exp (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	double fact = 1.0;
	for (k=0; k<1000; k++)
	{
		double term = tw / fact;
		acc += term* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		fact *= k+1;
		if (1.0e-14 > term) break;
	}

	return acc;
}

main (int argc, char *argv[])
{
	int i;

	// int nmax = 512;
	// int nmax = 431;
	int nmax = 19;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	for (i=0; i<nmax; i++)
	{
		double ts = 1.0;
		double x = i/((double)nmax);
		// double tw = takagi (w, x);
		// double ts = sin_takagi (w, x);
		// double tw = dtakagi (w, x);
		// double tw = log (takagi(w,x));
		// double tw = ttakagi(w,x);
		// double tw = div_takagi (w, x);
		// double tw = takagi_bumps (w, x);
		// double tw = takagi_prime (w, x);
		double tw = takagi_exp (w, x);


		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
	}
}
