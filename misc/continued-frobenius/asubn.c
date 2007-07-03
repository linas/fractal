/* 
 * asubn.c
 *
 * Graph of a_sub_n -- for manuscript
 * but also explore the general case better
 *
 * Linas Vepstas 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ache.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
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
	for (i=2; i<40; i++)
	{
// #define TRADITIONAL_A_SUB_N
#ifdef TRADITIONAL_A_SUB_N
		double x = a_sub_n (i);
		double r = sqrt (i*M_PI);
		r = exp (-2.0*r);
		r *= pow (i, -0.75);
		x /= r;
		double y = sqrt (4.0*M_PI*i) + 0.375*M_PI;
		y = -cos(y);
		y *= pow (0.5*M_PI, -0.25);
		printf ("%d	%8.6g	%8.6g\n", i, x, y);
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
		x *= exp (sqrt(4.0*M_PI*i));
		x *= pow (2.0/(M_PI*i*i*i), -0.25);

		double y = hurwitz_a_sub_n (i, m, k);
		y *= exp (sqrt(4.0*M_PI*i/((double) k)));
		y *= pow (2.0/(M_PI*i*i*i*k), -0.25);

		x = y;

		y = M_PI *((2.0*m/((double)k)) + 0.125);
		y -= sqrt (4.0*M_PI*i/((double) k));
		y = -sin (y);

		double z = 1.0/(y-prev);
		// z /= i*(i+1);
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
		prev = y;
#endif

#ifdef ETA_A_SUB_N
		double x = a_sub_n (i);
		x *= exp (sqrt (4.0*M_PI*i));
		double y = eta_a_sub_n (i);
		y *= exp (sqrt (2.0*M_PI*i));
		double z = gsl_sf_hzeta (0, ((double) m) / ((double)k));
x=z;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
#endif

// #define CHECK_B_SUB_N
#ifdef CHECK_B_SUB_N
		double y = b_sub_n (i);
		double z = b_sub_n_direct (i);
		printf ("%d	%8.6g   %8.6g	%8.6g\n", i, y,z, y-z);
#endif

#define RIEMANN_B_SUB_N
#ifdef RIEMANN_B_SUB_N
		double y = b_sub_n (i);
		y *= exp (sqrt(4*M_PI*i));
		// double z = -cos (sqrt(4*M_PI*i)+0.375*M_PI);
		double z = cos (sqrt(4*M_PI*i)-0.625*M_PI);
		z *= sqrt(sqrt(2*i/M_PI));
		// z = hurwitz_b_sub_n (i,1,1);
		// z *= exp (sqrt(4*M_PI*i));
		printf ("%d	%8.6g   %8.6g\n", i, y,z);
#endif

// #define HURWITZ_B_SUB_N
#ifdef HURWITZ_B_SUB_N
		double y = hurwitz_b_sub_n (i, m, k);
		double z = k*small_b_sub_n (i,m,k);
		printf ("%d	%8.6g   %8.6g\n", i, y,z);
#endif

// #define HURL_B_SUB_N
#ifdef HURL_B_SUB_N
		double y = hurwitz_b_sub_n (k, m, i);
		double z = i*small_b_sub_n (k,m, i);
		printf ("%d	%8.6g   %8.6g\n", i, y,z);
#endif

// #define WITZ_B_SUB_N
#ifdef WITZ_B_SUB_N
		double y = hurwitz_b_sub_n (i, 1, 1);
		y -= hurwitz_b_sub_n (i, 1, 2);
		y -= hurwitz_b_sub_n (i, 2, 2);
		double z = small_b_sub_n (i,1,1);
		z -= small_b_sub_n (i,1,2);
		z -= small_b_sub_n (i,2,2);
		printf ("%d	%8.6g   %8.6g\n", i, y,z);
#endif

// #define MUL_WITZ_B_SUB_N
#ifdef MUL_WITZ_B_SUB_N
		double y = hurwitz_b_sub_n (i, 1, 1);
		for (m=1; m<=k; m++)
		{
			y -= hurwitz_b_sub_n (i, m, k);
		}
		double z = small_b_sub_n (i,1,1);
		for (m=1; m<=k; m++)
		{
			z -= small_b_sub_n (i,m,k);
		}
		printf ("%d	%8.6g   %8.6g\n", i, y,z);
#endif

// #define DIGAMMA_CHECK
#ifdef DIGAMMA_CHECK
		double y = 0.0;
		k = i;
		for (m=1; m<=k; m++)
		{
			y += gsl_sf_psi (((double)m)/((double) k));
		}
		y += k *(M_EULER + log(k));
		printf ("%d	%8.6g\n", i, y);
#endif
	}
	return 0;
}
