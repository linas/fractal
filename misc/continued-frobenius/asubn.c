/* 
 * asubn.c
 *
 * Graph of a_sub_n -- for manuscript
 * but also explore the gneral case better
 *
 * Linas Vepstas 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ache.h"
#include <gsl/gsl_sf_zeta.h>

// ==========================================================
// return the harmonic numbers

long double harm_n2p1 (int n)
{
	int k;
	long double sum = 0.0L;
	for (k=1; k<=n; k++)
	{
		sum += 1.0L/((long double) k*(k+1));
	}
	return sum;
}

// ==========================================================
int
main (int argc, char * argv[])
{
	int i;

	if (3> argc)
	{
		fprintf (stderr, "Usage: %s <m> <k>\n", argv[0]);
		exit (1);
	}
	int m = atoi (argv[1]);
	int k = atof (argv[2]);

	double prev = 0.0;
	for (i=1; i<40; i++)
	{
#if TRADITIONAL_A_SUB_N
		double x = a_sub_n (i);
		double r = sqrt (i+1);
		r = exp (-4.0*r);
		x /= r;
		printf ("%d	%8.6g\n", i, x);
#endif
#if 0
		double x = t_sub_n (i, a);
		double y = a_sub_n (i);
		x -= 1.0/(2.0*a*((double)(i+1)));
		printf ("%d	%8.6g	%8.6g\n", i, x, y);
#endif
#ifdef HURWITZ_A_SUB_N
		double x = a_sub_n (i);
		double y = hz_a_sub_n (i, a);
		double z = 1.0/(y-prev);
		//z /= i*(i+1);
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
		prev = y;
#endif

// #define LFUNC_A_SUB_N
#ifdef LFUNC_A_SUB_N
		double x = a_sub_n (i);
		// double y = lfunc_a_sub_n (i, m, k);
		double y = eta_a_sub_n (i);
		// y += 2.0/((double)k) * harm_n2p1 (i);
		double z = 1.0/(y-prev);
		// z /= i*(i+1);
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
		prev = y;
#endif

#define ETA_A_SUB_N
#ifdef ETA_A_SUB_N
		double x = a_sub_n (i);
		x *= exp (sqrt (4.0*M_PI*i));
		double y = eta_a_sub_n (i);
		y *= exp (sqrt (2.0*M_PI*i));
		double z = gsl_sf_hzeta (0, ((double) m) / ((double)k));
x=z;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
#endif

#ifdef XXLFUNC_A_SUB_N
		double y = lfunc_a_sub_n (30, 1, i);
		printf ("%d	%8.6g\n", i, y);
#endif
	}
	return 0;
}
