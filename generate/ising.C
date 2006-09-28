
/*
 * ising.C
 *
 * Probability density of a two-sided 1D lattice model 
 * mapped to 1x1 square. Invariant under the Baker's 
 * transform.
 * 
 * One-dimensional Ising model (and Kac model too)
 * A given state of the system is encoded as 
 * a 2-adic string.  Encoding and general treatment
 *
 * Linas September 2005
 * Linas Vepstas Sept 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
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
	
	return 0.3 * s0*s1;
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
	return -1.0*en;
}

double qtent (double s)
{
	return tent(fquestion_mark (s));
}

// The farey/isola map
inline double pointy (double x)
{
	// double t = x - floor(x);
	double t = x;
	if (0.5 < t) return (1.0-t)/t;
	return t/(1.0-t);
}

/* analytic solution potential */
double vq (double s)
{
	if (0.5 < s) s  =1.0-s;
	
	double val = question_inverse (2.0*s);
	val = 2.0* log (1.0+val);
	val = log (2.0)-val;
	val = -val;
	return val;
}


/* Kac Model (which has shape of tent or cantor polynomial.) */
double kac (double s)
{
	// double lambda = 0.6666;
	double lambda = 0.5;
	// double lambda = 0.306852819;
	
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
	return (1.0-lambda)*s0*acc;
}

/* =========================================================== */
static double phase;

/* Return the finite-state energy of string s (length n) */
double energy (double (*interaction)(double), double s, int n)
{
	int i;

	double en = 0.0;
	for (i=0; i<n; i++)
	{
		double lam = cos (i*phase);
		en += lam * interaction (s);

		/* Shift one bit */
		if (s>= 0.5) s -= 0.5;
		s *= 2.0;
	}

	return en;
}

double cross_energy (double (*interaction)(double), double sl, double sr, int n)
{
	int i;

	double en = 0.0;
	for (i=0; i<n; i++)
	{
		/* Shift one bit */
		double s = 0.5*sr;
		sr = s;
		sl *= 2.0;
	  	if (1.0 < sl)
		{
			s += 0.5;
			sl -= 1.0;
		}	  
		double lam = cos ((i+1)*phase);
		en += lam * interaction (s);
		
	}

	return en;
}

/* Compute finite state partition */
double partition (double (*interaction)(double), double x, double y)
{
	int n = 10;
		
	// double y = fquestion_mark (x);
	// double en = energy (interaction, y, n);
	double en = 0.0;
	
	en = energy (interaction, x, n);
	en += cross_energy (interaction, y, x, n+1);
	
	en = exp (-en);
	return en;
}

/* explicitly for ising model only, to double-check that lattice
 * calcs are correct. Returns same answer as ising_density_alt. */
static double 
ising_density (double x, double y)
{
	int n = 10;

	double en = energy (nearest_neighbor, x, n);

	if (x<0.5 && y<0.5) en += 0.3;
	if (x<0.5 && y>=0.5) en -= 0.3;
	if (x>=0.5 && y<0.5) en -= 0.3;
	if (x>=0.5 && y>=0.5) en += 0.3;

	en += energy (nearest_neighbor, y, n);

	en = exp (-en);
	return en;
}

static double 
ising_density_alt (double x, double y)
{
	int n = 10;

	double en = energy (nearest_neighbor, x, n);
	y *= 0.5;
	if (x>0.5) y+= 0.5;
	en += energy (nearest_neighbor, y, n+1);

	en = exp (-en);
	return en;
}

static double 
density (double x, double y, int itermax, double param)
{
	phase = 0.25*M_PI*param;
	
	double p;
	// partition (pabola, n);
	// partition (vq, n);
	// partition (qtent, n);
	// partition (pointy, n);
 
	// p = partition (tent, x,y);
	// p = partition (kac, x,y);

	// should be equivalent to ising_density --- and it is.
	p = partition (nearest_neighbor, x,y);
	
	// p = ising_density (x,y);
	// p = ising_density_alt (x,y);
	return p;
}

DECL_MAKE_HEIGHT(density);

/* --------------------------- END OF LIFE ------------------------- */
