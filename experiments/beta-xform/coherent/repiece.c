/*
 * repiece.c
 * Failed attempt to build eigenfunctions as piecewise coherent states.
 * The code here is a hacked mess, because the basic idea can't work.
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

double gee(double beta, double omega, double y)
{
#define CONST
#ifdef CONST
	// const
	return 1.0;
#endif

// #define LINEAR
#ifdef LINEAR
	double off = -0.5 / (1.0 + omega*beta*beta);
	return y+off;
#endif

#ifdef QUAD
	double bee = -1.0 / (1.0 + omega * pow(beta, 3));
	double cee = -(0.25 + 0.5*bee) / (1.0 + omega * pow(beta, 2));
	return y*y + bee*y + cee;
#endif
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

double n2_psi_3(double beta, double omega, double alpha, double y)
{
	double m1 = 0.5*beta*(beta-1.0);
	if (m1 < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double n2_psi_2(double beta, double omega, double alpha, double y)
{
	if (0.5 < y) return 0.0;
	double m1 = 0.5*beta*(beta-1.0);
	if (y < m1) return 0.0;
	return psi(beta, omega, alpha, y);
}

double n2_psi_1(double beta, double omega, double alpha, double y)
{
	if (y < 0.5) return 0.0;
	if (0.5*beta < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double n2_psi(double beta, double omega, double y)
{
	if (0.5*beta < y) return 0.0;
	double a2 = beta * beta * omega;
	double a3 = beta * omega;
	double f1 = n2_psi_1(beta, omega, 1.0, y);
	double f2 = a2* n2_psi_2(beta, omega, beta, y);
	double f3 = a3* n2_psi_3(beta, omega, beta, y);
	double sum = f1 + f2 + f3;
	return sum;
}

// -----------------------------------------------------------

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s n k\n", argv[0]);
		exit (1);
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);

	// Poly roots
	double r1 = 0.5*(sqrt(5.0) + 1.0);
	double r2 = 1.465571231876768;

	double beta = 1.0;
	if (1 == n)
		beta = r1;
	if (2 == n)
		beta = r2;

	double omega = -pow(beta, -(k+2));
	if (2 == n)
		omega = -pow(beta, -(k+3));

	printf("#\n# beta=%f omega=%g\n#\n", beta, omega);

#define NPTS 2019
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
		double all = n2_psi(beta, omega, x);
		double p1 = n2_psi(beta, omega, x/beta);
		if (0.5*beta < x) p1=0.0;
		double p2 = n2_psi(beta, omega, x/beta + 0.5);
		double ell = (p1+p2) / beta;
		printf("%d	%f	%f	%f\n", j, x, all, ell);
		fflush(stdout);
	}
}
