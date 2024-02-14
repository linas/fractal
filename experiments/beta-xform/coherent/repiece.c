/*
 * repiece.c
 * Piecewise coherent states.
 *
 * February 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double modper(double period, double y)
{
	y -= period * floor (y / period);
	return y;
}

double ach(double y)
{
	return sin(44.0*y);
}

double gee(double beta, double omega, double y)
{
// #define CONST
#ifdef CONST
	// const
	return 1.0;
#endif

// #define LINEAR
#ifdef LINEAR
	double off = -0.5 / (1.0 + omega*beta*beta);
	return y+off;
#endif

// #define QUAD
#ifdef QUAD
	double bee = -1.0 / (1.0 + omega * pow(beta, 3));
	double cee = -(0.25 + 0.5*bee) / (1.0 + omega * pow(beta, 2));
	return y*y + bee*y + cee;
#endif

	if (y < 0.5) return ach(y);
	return -omega*beta*beta * ach(beta * (y-0.5));
}

// gee coherent function.
double psi(double beta, double omega, double alpha, double y)
{
	double m0 = 0.5*beta;
	if (m0 < y) return 0.0;

	double wn = 1.0;
	double yn = y;
	double sum = 0.0;
	while (1.0e-15 < fabs(wn))
	{
		double term = alpha * yn;
		sum += wn * gee(beta, omega, modper(m0, term));
		yn *= beta;
		if (m0 < yn) yn -= m0;
		wn *= omega;
	}
	return sum;
}

// -----------------------------------------------------------

double n1_psi_one(double beta, double omega, double alpha, double y)
{
	if (y < 0.5) return 0.0;
	if (0.5*beta < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double n1_psi_two(double beta, double omega, double alpha, double y)
{
	if (0.5 < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double n1_psi(double beta, double omega, double y)
{
	if (0.5*beta < y) return 0.0;
	double a2 = beta * omega;
	double f1 = n1_psi_one(beta, omega, 1.0, y);
	double f2 = a2* n1_psi_two(beta, omega, beta, y);
	double sum = f1 + f2;
	return sum;
}

// -----------------------------------------------------------

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s omega\n", argv[0]);
		exit (1);
	}
	double omega = atof(argv[1]);

	// Poly roots
	double r1 = 0.5*(sqrt(5.0) + 1.0);
	double beta = r1;
	// double omega = -pow(beta, -2);
	// double omega = -pow(beta, -4);

	printf("#\n# beta=%f omega=%g\n#\n", beta, omega);

#define NPTS 2019
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
		double all = n1_psi(beta, omega, x);
		double p1 = n1_psi(beta, omega, x/beta);
		if (0.5*beta < x) p1=0.0;
		double p2 = n1_psi(beta, omega, x/beta + 0.5);
		double ell = (p1+p2) / beta;
		printf("%d	%f	%f	%f\n", j, x, all, ell);
		fflush(stdout);
	}
}
