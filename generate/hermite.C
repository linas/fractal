
/*
 * hermite.C
 *
 * Green's function for the quantum harmonic oscillator
 * K(x,y)= sum_n H_n(x) H_N(y) / w_n   where w_n = n+1/2
 *
 * Linas Vepstas March 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/* compute hermite polynomial (physicists norm) up to order
 * nmax, using simple recursive algo, and place results in array "vals"
 */
void hermite_array (double x, int nmax, double *H_n)
{
	H_n[0] = 0.0;
	H_n[1] = 2.0*x;

	int n;
	for (n=1; n<nmax-1; n++)
	{
		H_n[n+1] = 2.0*(x*H_n[n] - n*H_n[n-1]);
	}
}

double **hermite_grid;
int nsteps;
int nmax;
double xmax;
double xdelta;
 
void init_hermite_grid (double exmax, double exdelta, int enmax)
{
	nsteps = (int) floor(xmax / xdelta) +1;
	nmax = enmax;
	xmax = exmax;
	xdelta = exdelta;

	hermite_grid = (double **) malloc (nsteps * sizeof (double *));
	double *arr = (double *) malloc (nsteps *nmax*sizeof (double));

	double x = 0.0;
	int i;
	for (i=0; i<nsteps; i++)
	{
		hermite_grid[i] = arr;
		hermite_array (x, nmax, arr);
		arr += nmax;
		x += xdelta;
	}
}

double hermite (double x, int n)
{
	double sgn = 1.0;
	if (x < 0.0)
	{
		x = -x;
		if (n%2) sgn = -1.0;
	}
	int iarr = (int) ((x+0.5*xdelta)/xdelta);
	return  sgn*hermite_grid[iarr][n];
}

static double hermite_green (double x, double y)
{
	double acc = 0.0;
	int n;
	double en = 0.5;
	double tn = 1.0;
	double nfac = 1.0;
	for (n=0; n<nmax; n++)
	{
		double term =  hermite (x,n) * hermite (y,n);
		term /= tn * nfac * en;
		acc += term;

		nfac *= n+1;
		en = n+0.5;
		tn *= 2.0;
		printf ("duude n=%d term=%g\n", n, term);
	}
	printf ("\n");

	acc *= exp (-0.5*(x*x+y*y));
	acc /= sqrt (M_PI);
	return acc;
}

static double hermite_green_wrap (double x, double y, int itermax, double param)
{
	int is_init = 0;
	if (!is_init)
	{
		is_init = 1;
		init_hermite_grid (9.0, 0.01 ,100);
	}
	return hermite_green (x,y);
}

DECL_MAKE_HISTO(hermite_green_wrap);

/* --------------------------- END OF LIFE ------------------------- */
