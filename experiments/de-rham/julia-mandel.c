/*
 * julia-mandel.c
 * Julia set of Mandelbrot iterated eqn.
 *
 * December 2023
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#define COMPLEX double complex

COMPLEX d0(COMPLEX zee, COMPLEX cee)
{
	return csqrt(zee-cee);
}

COMPLEX d1(COMPLEX zee, COMPLEX cee)
{
	return -csqrt(zee-cee);
}

COMPLEX iterate(double ex, COMPLEX cee)
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
		printf("%d	%g	%g\n", i, creal(zee), cimag(zee));
	}
}
