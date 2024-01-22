/*
 * unref.c
 *
 * Test patterns for unwrapped recursion relations for generalized
 * stretch-cut-stack map.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Invariant measure for beta shift (not beta transform)
double invar(double beta, double x)
{
	double midpnt = 0.5*beta;
	double obn = 1.0;
	double sum = 0.0;
	double norm = 0.0;
	for (int i=0; i<1000; i++)
	{
		if (x < midpnt) sum += obn;
		norm += midpnt*obn;

		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		if (obn < 1e-15) break;
	}
	return sum / norm;
}

double gp_invar(double beta, double x)
{
	return 0.5*beta*invar(beta, 0.5*beta*x);
}

// Eigenfunction for n=1 polynomial
// Use with beta= 1.618034
double n1(double beta, double x)
{
	double whack = 1.105;  // some magic normalization!?
	double m0 = 0.5*beta;
	if (x < 0.5) return whack*(beta*x-0.5);
	if (x < m0) return whack* (x-0.5);
	return 0.0;
}

double gp_n1(double beta, double x)
{
	return 0.5*beta*n1(beta, 0.5*beta*x);
}

// Eigenfunction for n=1 polynomial w/ quadratic generator
double quad_n1(double beta, double x)
{
	double whack = 1.0;  // some magic normalization!?
	double m0 = 0.5*beta;
	if (x < 0.5) return whack*(beta*x*x-x+0.125);
	if (x < m0) return whack*(x*x-x+ beta*0.125);
	return 0.0;
}

double gp_quad_n1(double beta, double x)
{
	return 0.5*beta*quad_n1(beta, 0.5*beta*x);
}

// Eigenfunction for n=2 polynomial
// Use with beta= 1.4655712318
double n2(double beta, double x)
{
	double whack = 1.13;  // some magic normalization!?
	double m0 = 0.5*beta;
	double m1 = 0.5*beta * (beta-1.0);
	if (x < m1) return whack*(beta*beta*x - 0.5);
	if (x < 0.5) return whack*(beta*x-0.5);
	if (x < m0) return whack*(x-0.5);
	return 0.0;
}

double gp_n2(double beta, double x)
{
	return 0.5*beta*n2(beta, 0.5*beta*x);
}

// Eigenfunction for n=3 polynomial
// Use with beta= 1.839286755
double n3(double beta, double x)
{
	double m0 = 0.5*beta;
	double m1 = 0.5*beta * (beta-1.0);
	double whack = 0.95;  // some magic normalization!?
	if (x < 0.5) return whack*m1*(beta*x- m1 + 0.25/beta);
	if (x < m1) return whack*m1*((beta+1.0)*x/beta - m1);
	if (x < m0) return whack*m1*(x- m1 + 0.25/beta);
	return 0.0;
}

double gp_n3(double beta, double x)
{
	return 0.5*beta*n3(beta, 0.5*beta*x);
}

// ==============================================================
