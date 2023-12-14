/*
 * minside.c
 *
 * Explore rate of convergence inside the mandelbrot set.
 *
 * December 2023
 */

#include "brat.h"

#include <iostream>
#include <complex.h>
#define COMPLEX std::complex<double>
// #define COMPLEX double complex

static double mandelbrot_convergence(double re_q, double im_q, int itermax, double param)
{
	COMPLEX cee = re_q + I * im_q;
	COMPLEX zee = cee;
	COMPLEX *orb = (COMPLEX *) malloc (itermax * sizeof(COMPLEX));

#define EPS 1e-2
	for (int i=0; i<itermax; i++)
	{
		orb[i] = zee;

		zee = zee*zee+cee;
		if (1/EPS < abs(zee))
			return i;
	}

	int cnt = 0;
	return cnt;
}

DECL_MAKE_HEIGHT(mandelbrot_convergence);

/* --------------------------- END OF LIFE ------------------------- */
