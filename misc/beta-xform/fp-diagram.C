/*
 * fp-diagram.C
 *
 * Feigenbaum-style diagram for the undermap, but computed
 * "exactly", not statisically, so free of noise. I'm calling
 * this the Frobenius-Perron diagram.
 *
 *
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bitops.h"
#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 * Non-renomralized 
 * Recursve frobenius-perron eigenfunction for the undershift
 */
double reig(double x, double K, int niter)
{
	double tkay = 2.0*K;

	if (K < x) return 0.0;
	if (niter < 0) return 1.0;
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, niter-1) / tkay;
	}
	double sum = reig(x/tkay, K, niter-1);
	sum += reig(0.5 + x/tkay, K, niter-1);
	sum /= tkay;
	return sum;
}

/*
 * Normalization for 
 * Recursive frobenius-perron eigenfunction for the undershift
 * That is, the integral of the unrenomralized function.
 */
double fpeig(double K, int niter)
{
#define NPTS 2803
	double acc = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NPTS);
		double y = reig(x, K, niter);
		acc += y;
	}
	acc /= NPTS;
	return acc;
}

double niter = 18;

static void bifurcation_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double K,
                                 int itermax,
                                 double omega)
{
#if 0
	double norm = fpeig (K, itermax);

	for (int j=0; j<array_size; j++)
	{
		double x = (((double) j) + 0.5) / ((double) array_size);
		double eig = reig(x, K, itermax);
		array[j] = eig / norm;
	}
#endif

#if 0
if (K < 0.9) itermax *= 1.2;
if (K < 0.8) itermax *= 1.6;
if (K < 0.7) itermax *= 1.5;
if (K < 0.65) itermax *= 1.4;
if (K < 0.6) itermax *= 1.8;
if (K < 0.55) itermax *= 2;
#endif

int perow = itermax;
double bump = 1.05;
if (K < 0.6) { perow *= 2; bump = 1.1; }
if (K < 0.55) { perow *= 4; bump = 1.2; }

itermax = niter;

time_t start = time(0);
	for (int j=0; j<array_size; j++)
	{
		double x = (((double) j) + 0.1) / ((double) array_size);
		double eig = reig(x, K, itermax);
		array[j] = eig;

		x = (((double) j) + 0.3) / ((double) array_size);
		eig = reig(x, K, itermax);
		array[j] += eig;

		x = (((double) j) + 0.5) / ((double) array_size);
		eig = reig(x, K, itermax);
		array[j] += eig;

		x = (((double) j) + 0.7) / ((double) array_size);
		eig = reig(x, K, itermax);
		array[j] += eig;

		x = (((double) j) + 0.9) / ((double) array_size);
		eig = reig(x, K, itermax);
		array[j] += eig;
	}
time_t end = time(0);
printf("end %g itermax=%d niter=%f time=%lu\n", K, itermax, niter, end-start);
if (end-start < perow/2) niter *= bump;
if (end-start < perow) niter *= bump;
if (4*perow < end-start) niter /= 1.02;


	double norm = 0.0;
	for (int j=0; j<array_size; j++)
	{
		norm += array[j];
	}
	norm /= array_size;

	for (int j=0; j<array_size; j++)
	{
		array[j] /= norm;
	}
}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
