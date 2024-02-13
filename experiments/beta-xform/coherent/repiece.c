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
	return 1.0;
	double off = -0.5 / (1.0 + omega*beta*beta);
	return y+off;
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

double psi_one(double beta, double omega, double alpha, double y)
{
	if (y < 0.5) return 0.0;
	if (0.5*beta < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double psi_two(double beta, double omega, double alpha, double y)
{
	if (0.5 < y) return 0.0;
	return psi(beta, omega, alpha, y);
}

double psi_all(double beta, double omega, double y)
{
	if (0.5*beta < y) return 0.0;
	double a2 = beta * omega;
	double f1 = psi_one(beta, omega, 1.0, y);
	double f2 = a2* psi_two(beta, omega, beta, y);
	double sum = f1 + f2;
	return sum;
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s omega k\n", argv[0]);
		exit (1);
	}
	double omega = atof(argv[1]);
	int k = atoi(argv[2]);

	double beta = 0.5*(sqrt(5.0) + 1.0);

	printf("#\n# beta=%f k=%d omega=%g\n#\n", beta, k, omega);

#define NPTS 2019
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
		double all = psi_all(beta, omega, x);
		double p1 = psi_all(beta, omega, x/beta);
		if (0.5*beta < x) p1=0.0;
		double p2 = psi_all(beta, omega, x/beta + 0.5);
		double ell = (p1+p2) / beta;
		printf("%d	%f	%f	%f\n", j, x, all, ell);
		fflush(stdout);
	}
}
