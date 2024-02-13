/*
 * repiece.c
 * Build eigenfunctions as piecewise coherent states.
 *
 * February 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Sawtooth for linear g.
double sgen(double beta, double alpha, double y)
{
	double m0 = 0.5*beta;
	double gen = alpha*y;
	gen -= m0 * floor (gen / m0);
	return gen;
}

// Line-generated coherent function.
double coh(double beta, double omega, double alpha, double y)
{
	double m0 = 0.5*beta;
	if (m0 < y) return 0.0;

	double wn = 1.0;
	double yn = y;
	double sum = 0.0;
	while (1.0e-15 < wn)
	{
		double term = alpha * yn;
		sum += wn * (term - m0*floor(term/m0)); // sgen()
		yn *= beta;
		if (m0 < yn) yn -= m0;
		wn *= omega;
	}
	return sum;
}

double goldcoh(double omega, int k, double y)
{
	double beta = 0.5*(sqrt(5.0) + 1.0);
	double alpha = beta * k;

	// Interval 2
	if (y < 0.5)
		return coh(beta, omega, alpha, y);

	// Interval 1
	if (y < 0.5*beta)
		return -coh(beta, omega, alpha, y/beta) / (omega*beta*beta);

	return 0.0;
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

#if LATER
#define NPTS 2019
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
		double y = goldcoh(omega, k, x);
		double bra1 = goldcoh(omega, k, x/beta);
		if (0.5*beta < x) bra1 = 0.0;
		double bra2 = goldcoh(omega, k, x/beta+0.5);
		if (0.5 < x) bra2 = 0.0;
		printf("%d	%f	%f	%f	%f\n", j, x, y, bra1, bra2);
		fflush(stdout);
	}
#endif

#define SANITY_CHECK
#ifdef SANITY_CHECK
	// Quadruple-check the coherent sum. Yes, it behaves exactly how
	// it should. So this is not where the problem lies.
	double alpha = k*beta;
#define NPTS 19
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
#ifdef LCHECK
		double y = coh(beta, omega, alpha, x);
		double bra1 = coh(beta, omega, alpha, x/beta);
		double sha1 = k*x + omega * y - bra1;
		printf("%d	%f	%f	%f	%g\n", j, x, y, bra1, sha1);
#endif
// #define RCHECK
#ifdef RCHECK
		double y = coh(beta, omega, alpha, x);
		double bra2 = coh(beta, omega, alpha, x/beta + 0.5);
		double gen = sgen(beta, alpha, x/beta + 0.5);
		double sha2 = gen + omega * y - bra2;
		printf("%d	%f	%f	%f	%g\n", j, x, y, bra2, sha2);
#endif

#ifdef TERM_COLLECT
		// Passes great
		double a2 = -1.0 / beta;
		double ss = a2*sgen(beta, alpha, x/beta);
		ss += sgen(beta, alpha, x/beta + 0.5);
		ss += a2 * sgen(beta, alpha, x/beta) / beta;
		printf("%d	%f	%g\n", j, x, ss);
#endif

#ifdef TERM_CONDENSE
		// Passes great
		double a2 = -1.0 / beta;
		double ss = a2*beta*sgen(beta, alpha, x/beta);
		ss += sgen(beta, alpha, x/beta + 0.5);
		printf("%d	%f	%g\n", j, x, ss);
#endif

#ifdef LAST_TERM
		// Passes great
		double a2 = -1.0 / beta;
		double ss = a2*beta*sgen(beta, alpha, x/beta);
		ss += sgen(beta, alpha, x/beta);
		printf("wat	%d  %f	%g\n", j, x, ss);
#endif

#ifdef TOP_ROW
		// works
		double y = goldcoh(omega, k, x);
		double top = omega*beta*y + goldcoh(omega, k, x/beta)/beta;
		printf("%d	%f	%f	%g\n", j, x, y, top);
#endif

#ifdef BOT_ROW
// borken
		double y = goldcoh(omega, k, x);
		double a2 = -1.0 / beta;
		double bot = -omega*y;
		bot -= -goldcoh(omega, k, x/beta) / beta;
		bot -= goldcoh(omega, k, x/beta + 0.5);
		printf("%d	%f	%f	%g\n", j, x, y, bot);
#endif

		fflush(stdout);
	}
#endif
}
