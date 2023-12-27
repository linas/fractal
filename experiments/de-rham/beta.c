/*
 * beta.c
 * Crazy experiments.
 *
 * December 2023
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define COMPLEX double complex

#if 1
COMPLEX d0(COMPLEX zee, COMPLEX cee)
{
	return creal(cee) * zee;
}

COMPLEX d1(COMPLEX zee, COMPLEX cee)
{
	return cimag(cee) * zee / (1+zee);
}
#endif

// This iteration is backwards from the deRham defintion...
// It reverses the order of function composition.
COMPLEX reverse_iterate(double ex, COMPLEX cee)
{
	COMPLEX zee = 0.0;
	for (int i=0; i<30; i++)
	{
		if (ex < 0.5)
			zee = d0(zee, cee);
		else
		{
			ex -= 0.5;
			zee = d1(zee, cee);
		}
		ex *= 2.0;
	}

	return zee;
}

// This iteration gives a continuous curve, by iterating the
// least-significant bits first.
COMPLEX iterate(double ex, COMPLEX cee)
{
	int digits = 30;
	double bigx = ex * (1UL << digits);
	long nx = floor(bigx);

	COMPLEX zee = 1.0;
	for (int i=0; i<digits; i++)
	{
		if (nx & 1UL)
			zee = d1(zee, cee);
		else
			zee = d0(zee, cee);

		nx >>= 1;
	}
	return zee;
}

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <re> <im> <npts>\n", argv[0]);
		exit(1);
	}
	double re = atof(argv[1]);
	double im = atof(argv[2]);
	int npts = atoi(argv[3]);
	COMPLEX cee = re + I*im;
	for (int i=0; i< npts; i++)
	{
		double x = (((double) i) + 0.5) / ((double) npts);
		COMPLEX zee = iterate(x, cee);
		printf("%d	%g	%g	%g\n", i, x, creal(zee), cimag(zee));
	}
}
