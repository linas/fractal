/*
 * moments.c
 *
 * Compute moments w.r.t. a measure that consists of Dirac deltas.
 * One delta per midpoint, midpoint k weighted with 1/beta^k.
 * This measue makes computation of moments and midpoints very
 * easy and fast.
 *
 * What did I learn from this? Nothing in particular. What can I do
 * with these polynomials? Nothing I can think of.
 *
 * January 2024
 */

#include <math.h> 
#include <stdio.h>
#include <stdlib.h> 

// Hausdorff moment (Hamburger moment on unit interval)
// This is an integral of the monmial x^n times a sequence
// of Dirac delta funcs, located at midpoints, and weighted
// with 1/beta^k for midpoint k. Thus, the integral becomes
// a sum, and a rapidly converging one, at that.
double fmoment(double beta, int n)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow(tk, n);
		sum += tkn / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

double coeff[20][20];

// Evaluate polynomial n at location x.
double orthonormo(int n, double x)
{
	double sum = 0.0;
	for (int j=0; j<n; j++)
		sum += coeff[n][j] * orthonormo (j, x);
	sum += coeff[n][n] * pow (x, n);
	return sum;
}

// Compute the integral of polynomial n times polynomial m
// Since these should be orthonormal, this should return
// zero or one, always. Since the measure is just a sequence
// of delta funcs, the integral is easy & fast to compute.
double prod(double beta, int n, int m)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double pa = orthonormo(n, tk);
		double pb = orthonormo(m, tk);
		sum += pa * pb / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

// Compute overlap between monomial x^n and polynomial p_j
// Used for Gaussian elimination & orthogonalization.
double ortho(double beta, int n, int j)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow(tk, n);
		double ort = orthonormo(j, tk);
		sum += ort * tkn / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

// Compute polynomial coefficients up to order max. Once
// this is done, the polynomials can be evaluated.
void setup(double beta, int max)
{
	double msq = fmoment(beta, 0);
	coeff[0][0] = 1.0 / sqrt(msq);

	for (int n=1; n< max; n++)
	{
		for (int j=0; j<n; j++)
			coeff[n][j] = -ortho(beta, n, j);
		coeff[n][n] = 1.0;
		double msq = prod(beta, n, n);
		double rms = 1.0 / sqrt(msq);
		for (int j=0; j<=n; j++)
			coeff[n][j] *= rms;
	}
}

// ==============================================================

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
#if 0
for (int n=0; n<nset; n++) {
for (int m=0; m<=n; m++) {
printf("%d %d co= %f  rat= %f\n", n, m, coeff[n][m], coeff[n][m]/coeff[n][n]);
}
printf("---\n");
}
#endif

	int nset = 8;
	setup(beta, nset);
	int imax = 432;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		printf("%d	%f", i, x);
		for (int n=0; n<nset; n++)
		{
			double poly = orthonormo(n, x);
			printf("	%f", poly);
		}
		printf("\n");
		fflush(stdout);
	}

#ifdef PRINT_MOMENTS
	int imax = 1200;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;

		printf("%d	%f", i, beta);
		for (int n=0; n<8; n++)
		{
			double fmom = fmoment(beta, n);
			printf("	%.10f", fmom);
		}
		printf("\n");
		fflush(stdout);
	}
#endif
}