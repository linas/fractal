
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
	nsteps = xmax / xdelta +1;
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
	int iarr = x/xdelta;
	return  hermite_grid[iarr][n];
}

static double hermite_green (double x, double y)
{
	double acc = 0.0;
	int n;
	for (n=0; n<nmax; n++)
	{
	}
	return acc;
}

static double hermite_green_wrap (double x, double y, int itermax, double param)
{
	return hermite_green (x,y);
}

DECL_MAKE_HISTO(hermite_green_wrap);

/* --------------------------- END OF LIFE ------------------------- */
