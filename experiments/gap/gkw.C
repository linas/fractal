
/* gkw.C
 *
 * Explore gauss-kuzmin-wirsing symmetries
 *
 * Linas Vepstas November 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

double sumq (double x)
{
	ContinuedFraction f;
	f.SetEvenize();

	double acc = 0.0;
	double tn = 0.5;
	int n;
	for (n=1; n<30; n++)
	{
		double arg = 1.0 / (((double) n) +x);
		arg = arg*arg*arg*arg;
		f.SetReal (arg);
		double term = f.ToFarey();
		term *= tn;
		acc += term;
		tn *= 0.5;
	}
	return acc;
}

double qpsi (double x)
{
	ContinuedFraction f;
	f.SetEvenize();

	double acc = 0.0;
	int n;
	for (n=1; n<1000; n++)
	{
		double arg = 1.0 / (((double) n) +x);
		arg = arg*arg*arg*arg;
		acc += arg;
	}
	f.SetReal (acc);
	double qp = f.ToFarey();

	return qp;
}

main (int argc, char *argv[])
{
	ContinuedFraction f;
	f.SetEvenize();
	int i;

	int nmax = 431;
	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/ ((double) nmax);

		double qp = qpsi (x);
		double sq = sumq (x);
		sq *= 2.0;
		
		f.SetReal (x);
		double q = f.ToFarey();
		q = 1.0-q;

		printf("%5d	%8.6g	%8.6g	%8.6g	%8.6g\n", i,x, qp, sq, q);

	}
}

