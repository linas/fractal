/*
 * kac.c
 *
 * One-dimensional Kac model. General resmblance to
 * the code in ising.c except that more careful attention 
 * is given to the cylinder-set topology.
 * 
 * A cylinder set is encoded as a 2-adic string (a real)
 * with the number of bits in the string specified.
 *
 * Linas September 2005
 * Linas September 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "question.h"

/* Interaction functions: Given a string 's' of length 'n'
 * encoding a given cylinder set, * return the interaction 
 * energy associated with that state. The string is represented 
 * as a real number between zero and one: it is the string of 
 * binary bits of the expansion of the real in binary. Only
 * the first 'n' bits are used.
 */

/* Nearest neighbor interaction */
double nearest_neighbor (double s, int len)
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

/* Kac Model (which has shape of tent or cantor polynomial.) */
double kac (double s, int len)
{
	// double lambda = 0.6666;
	double lambda = 0.5;

	if (0 >= len) return 0.0;
	
	double s0 = -1.0;
	if (s>= 0.5) {
		s0 = 1.0;
		s -= 0.5;
	}
	s *= 2.0;
	
	double lp = lambda;
	double acc = 0.0;
	int i;
	for (i=0; i<len; i++)
	{
		double s1 = -1.0;
		if (s>= 0.5) {
			s1 = 1.0;
			s -= 0.5;
		}
		s *= 2.0;
	
		acc += lp * s1;
		lp *= lambda;
	}
	return -(1.0-lambda)*s0*acc;
}

/* 
 * Return the energy associated with a cylinder set encoded
 * by a string s of length n.
 */
double energy (double (*interaction)(double, int), double s, int len)
{
	int i;

	double h_n = 0.0;
	for (i=0; i<len; i++)
	{
		h_n += interaction (s, len-i);

		/* Shift one bit */
		if (s>= 0.5) s -= 0.5;
		s *= 2.0;
	}

	return h_n;
}

/* Compute finite state partition */

double partition (double (*interaction)(double, int), int n)
{
	double ez = 0.0;
	double z = 0.0;

	int m = (1<<n);
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
