
/* 
 * ising.C
 *
 * Ising model inspired takagi curve stuff 
 *
 * Linas October 2004
 * Linas September 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
ContinuedFraction far;

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

/* The main, core basic takagi curve */
long double takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<50; k++)
	{
		// long double term = tw* balanced_triangle (tp*x);
		long double term = tw* triangle (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* x = p / 2**n */
int itriangle (int p, int n)
{
	unsigned int deno = 1<<n;

	if (p >= deno) p-= deno;
	if (p == 0) return 0;

	if (2*p < deno) return p;

	p = deno - p;
	return p;
}

/* x = p / 2**n */
long double itakagi (long double w, int p, int n)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<n; k++)
	{
		p = itriangle (p, n-k);
		if (0 == p) break;
		long double term = 2.0L * tw * ((long double) p) / ((long double) (1<<(n-k)));
		acc += term;
		tp *= 2.0L;
		tw *= w;
	}

	return acc;
}

// ===============================================
// ===============================================

int
main (int argc, char *argv[])
{
	ContinuedFraction fr;

	int i;

	if (argc <4)
	{
		printf ("Usage: %s <log2-len> <w-value> <exp-base>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	double w = atof(argv[2]);
	double base = atof(argv[3]);
	
	// int nmax = 512;
	// int nmax = 432;
	// int nmax = 1717;
	// int nmax = 2048;

	int nmax = 1<<p;
	
	double scale = sqrt ((double) nmax);
	printf ("#\n# scale=%g\n#\n", scale);

	double fa = 0.0;
	double accr = 0.0;
	double acci = 0.0;
	double xlast = 0.0;
	double step = 0.0012;

	int ia = 0;
#define ARRSZ 4000
	double xa[ARRSZ];
	double twa[ARRSZ];
	double accar[ARRSZ];
	double accai[ARRSZ];
	double faa[ARRSZ];
	double norm = 1.0;

	srand (p);

	double delta = 1.0 / ((double)nmax);
	double x = 0.0;
	for (i=1; i<nmax+30; i++)
	{
		// x = i* delta;

		/* jitter, as this can make a difference */
		double dh = ((double) rand()) / ((double)RAND_MAX);
		dh *= 2.0*delta;
		x += dh;
		
		// double ts = isola (w, x);
		double tw = takagi (w, x);
		// double tw = itakagi (w, i, p);
		
		// double ts = exp (-tw);
		// double ts = pow (4.0, -tw);
		// double ts = pow (M_PI, -tw);
		// double ts = pow (base, -tw);
		double ts = cos (base*tw);
		accr += ts* dh;
		ts = sin (base*tw);
		acci += ts* dh;

		if (x <= 0.5) norm = accr;
		if (x > 1.0) break;

		if (x> xlast + step) {
			fa = InvFarey(x);
			xlast += step;
			printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, accr, acci, fa);
			xa[ia] = x;
			twa[ia] = tw;
			accar[ia] = accr;
			accai[ia] = acci;
			faa[ia] = fa;
			ia ++;
		}

		// fflush (stdout);
	}
return 0;

	norm = 2.0*norm;
	norm = 1.0/norm;
	for (i=0; i<ia; i++)
	{
		accar[i] *= norm;
		accai[i] *= norm;
		// fr.SetReal (acca[i]);
		// double y = fr.ToFarey();
		// faa[i] = xa[i] - y;
		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, xa[i], twa[i], accar[i], accai[i], faa[i]);
	}
	printf ("#\n# last= %f\n", accr);
	fprintf (stderr, "#\n# last= %g\n", accr);
	fprintf (stderr, "#\n# norm= %g\n", norm);
}
