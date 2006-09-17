/*
 * ising.c
 *
 * 1 dimensional Ising model (and Kac model too)
 * A given state of the system is encoded as 
 * a 2-adic string.
 *
 * Linas September 2005
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "question.h"

/* Interaction functions: Given a string 's' encoding 
 * a given state, return the interaction energy associated 
 * with that state. The string is represented as a real
 * number between zero and one: it is the string of binary 
 * bits of the expansion of the real in binary.
 */

/* Nearest neighbor interaction */
double nearest_neighbor (double s)
{
	double s0,s1;

	s0 = -1.0;
	if (s>= 0.5) {
		s0 = 1.0;
		s -= 0.5;
	}
	s *= 2.0;
	
	s1 = -1.0;
	if (s>= 0.5) s1 = 1.0;
	
	return -0.6 * s0*s1;
}

double pabola (double s)
{
	return -0.2 * s * (1.0-s);
}

double tent (double s)
{
	// double en = (s>0.5)? 2.0*s : 2.0*(1.0-s);
	double en = (s>0.5)? 2.0*(1.0-s) : 2.0*s;
	en -= 0.5;
	return 1.0*en;
}

// The farey/isola map
inline double pointy (double x)
{
	// double t = x - floor(x);
	double t = x;
	if (0.5 < t) return (1.0-t)/t;
	return t/(1.0-t);
}


/* Kac Model (which has shape of tent or cantor polynomial.) */
double kac (double s)
{
	double lambda = 0.6666;
	// double lambda = 0.5;
	
	double s0 = -1.0;
	if (s>= 0.5) {
		s0 = 1.0;
		s -= 0.5;
	}
	s *= 2.0;
	
	double lp = lambda;
	double acc = 0.0;
	while (1)
	{
		double s1 = -1.0;
		if (s>= 0.5) {
			s1 = 1.0;
			s -= 0.5;
		}
		s *= 2.0;
	
		acc += lp * s1;
		lp *= lambda;

		if (lp < 1.0e-18) break;
	}
	return -(1.0-lambda)*s0*acc;
}

/* Return the finite-state energy of string s (length n) */
double energy (double (*interaction)(double), double s, int n)
{
	int i;

	double en = 0.0;
	for (i=0; i<n; i++)
	{
		en += interaction (s);

		/* Shift one bit */
		if (s>= 0.5) s -= 0.5;
		s *= 2.0;
	}

	return en;
}

/* Compute finite state partition */

double partition (double (*interaction)(double), int n)
{
	double ez = 0.0;
	double z = 0.0;

	int m = (1<<n) +43;
	int prt = m/24000;
	int i;

	double om = 1.0 / ((double) m);
	double delta = om;
	double sprev = 0.0;
	
	for (i=1; i<m; i++)
	{
		double x = om * ((double) i);
		double y = question_mark (i,m);
		
		double en = energy (interaction, y, n);
		
		en = exp (en);
		z += delta * en;

		if (i%prt == 0) {
			printf ("%d	%10.8g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, y, interaction(x), en, z);
			fflush (stdout);
		}
	}

	printf ("# partition = %g\n", z);

	return z;
}

int
main (int argc, char * argv[]) 
{
	int n = 10;

	if (argc < 2) {
		fprintf (stderr, "Usage: %s <n>\n", argv[0]);
		exit (1);
	}
	n = atoi (argv[1]);

	printf ("#\n# n=%d\n#\n",n);
	// partition (nearest_neighbor, n);
	// partition (pabola, n);
	 partition (tent, n);
	// partition (kac, n);
	// partition (pointy, n);
}
