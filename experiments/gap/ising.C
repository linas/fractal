
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
#include "question.h"
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
		// long double term = tw* fabs(sin (tp*x*M_PI));
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

// ===============================================
// A simple minded integer variation.  Performs
// only a finite sum
// 
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
// Try to perform the full summation, in pure integer format.
//  x = p / q;
// 

long double iexact_takagi (long double w, long long p, long long q)
{
	long long pee = p;

	if (0 == p)
	{
		fprintf (stderr, "P is zero!! \n");
		return 0.0L;
	}
	
	/* first,sum up the transient bit */
	long double tracc = 0.0L;
	long double tw = 1.0L;
	while(1)
	{
		p *= 2;
		if (p >= q) break;
		tracc += tw * ((long double) p);
		tw *= w;
	}

	/* make note of the start at the peak. */
	long long pstart = p/2;
	int len = 1;
	long double tran_tw = tw;
	
	/* one term by hand */
	p = 2*q-p;
	long double repacc = ((long double) p);
	tw = w;
	
	/* Now, sum up the repeating part */
	while (p != pstart)
	{
		p *= 2;
		if (p >= q) 
		{
			p = 2*q-p;
		}
		repacc += tw * ((long double) p);
		len ++;
		tw *= w;
		if (len > 60) 
		{
			fprintf (stderr, "oh no len too long!!!!\n");
			break;
		}
	}
printf ("x=%lld/%lld len=%d\n", pee,q, len);

	/* Result is transient plus repeating */
	long double acc = tracc + tran_tw * repacc / (1.0L - tw);

	acc /= (long double) q;
	return acc;
}

// ===============================================
// ===============================================

void
do_integral (int p, double w, double base)
{
	ContinuedFraction fr;
	int i;
	
	// int nmax = 432;
	int nmax = 1<<p;
	
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
	for (i=1; x < 1.0; i++)
	{
		// x = i* delta;

		/* jitter, as this can make a difference */
		double dh = ((double) rand()) / (((double)RAND_MAX)+1.0);
		dh *= 2.0*delta;
		x += dh;
		
		// double ts = isola (w, x);
		// double tw = takagi (w, x);
		double q = 3.0 + 458.0 * ((double) rand()) / (((double)RAND_MAX)+1.0);
		long long deno = (long long ) q;
		long long num = (long long) (x* ((double) deno));
		double tw = iexact_takagi (w, num, deno);
		// double tw = itakagi (w, i, p);
		
		// double ts = exp (-tw);
		// double ts = pow (4.0, -tw);
		// double ts = pow (M_PI, -tw);
		double ts = pow (base, -tw);
		accr += ts* dh;
#ifdef POLAR_PLOT
		double ts = cos (base*tw);
		accr += ts* dh;
		ts = sin (base*tw);
		acci += ts* dh;
#endif

		if (x <= 0.5) norm = accr;
		if (x > 1.0) break;

if(x>0.5) {
	printf ("%g	%g\n", w, norm);
	fflush (stdout);
	return;
}

		if (x> xlast + step) {
			fa = question_inverse(x);
			xlast += step;
			// printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, accr, acci, fa);
			xa[ia] = x;
			twa[ia] = ts;
			accar[ia] = accr;
			accai[ia] = acci;
			faa[ia] = fa;
			ia ++;
		}

		// fflush (stdout);
	}

	norm = 2.0*norm;
	norm = 1.0/norm;
	for (i=0; i<ia; i++)
	{
		accar[i] *= norm;
		accai[i] *= norm;
		// fr.SetReal (acca[i]);
		// double y = fr.ToFarey();
		// faa[i] = xa[i] - y;
		// printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, xa[i], twa[i], accar[i], accai[i], faa[i]);
	}
	printf ("#\n# last= %f\n", accr);
	fprintf (stderr, "#\n# last= %g\n", accr);
	fprintf (stderr, "#\n# norm= %g\n", norm);
}

int
main (int argc, char *argv[])
{
	if (argc <4)
	{
		printf ("Usage: %s <log2-len> <w-value> <exp-base>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	double w = atof(argv[2]);
	double base = atof(argv[3]);
	
	// do_integral (p,w,base);
	
	for (w=0.8; w < 1.0; w+= 0.005) { 
		do_integral (p,w,base);
	}
}
