/*
 * almost.c
 * Explore the almost-eigenfunctions that are coherent states
 * built from the Renyi-Parry bitstream
 *
 * Dec 2018
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// My downshift; beta=2K so T(x) = case bx or b(x-0.5)
double downshift(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

// Standard beta-xform is t(x) = (b x) mod 1
double beta_xform(double x, double beta)
{
	double prod = x * beta;
	return prod - floor(prod);
}

// Return the Renyi-Parry bit d_n (but for shift, not for xform.)
int d_n(int n, double x, double beta)
{
	double tn = 0.5*beta;
	if (tn < x) return 0;

	// compute T^N(b/2)
	for (int i=0; i<n; i++)
	{
		tn = downshift(tn, 0.5*beta);
	}
	return (x<tn);
}

// Build contorted invariant measure. i.e. stick a z^n into the sum.
double complex rhoz(double x, double beta, complex double z)
{
	if (0.5*beta < x) return 0.0;

	double ob = 1.0/beta;
	double obn = 1.0;
	complex double zn = 1.0;

	// accumulated sum
	complex double rho = 0.0;

	double tn = 0.5*beta;

	int cnt=0;
	while (1.0e-16 < obn*cabs(zn))
	{
		if (x < tn) {
			rho += obn * zn;
		}

		// compute 1/beta^N
		obn *= ob;

		// compute T^N(b/2)
		tn = downshift(tn, 0.5*beta);

		// compute z^n;
		zn *= z;

		cnt++;
		if (3000< cnt) break;
	}

	return rho;
}

// ================================================================

// Print the contorted density
void print_vee(double beta, double complex zee)
{

#define NPTS 3801
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		double complex rho = rhoz(x, beta, zee);
		double re = creal(rho);
		double im = cimag(rho);

		// beta shift applied to above
		double complex shift_rho = 0.0;
		if (x < 0.5*beta)
		{
			shift_rho = rhoz(x/beta, beta, zee);
			shift_rho += rhoz(x/beta + 0.5, beta, zee);
			shift_rho /= beta;
		}

		// rho divided by z;
		double complex rdz = rho / zee;

		// difference
		double complex cnst = shift_rho - rdz;
		double rec = creal(cnst);
		double imc = cimag(cnst);
		printf("%d	%g	%g	%g	%g	%g\n", i, x, re, im, rec, imc);
	}
}

// ================================================================
void print_d(int n, double beta)
{
#undef NPTS
#define NPTS 201
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		int dn = d_n(n, x, beta);
		int dnlo = d_n(n, x/beta, beta);
		int dnhi = d_n(n, x/beta+0.5, beta);
		printf("%d	%8.6f	%d	%d	%d	%d\n", i, x, n, dn, dnlo, dnhi);
	}
}

// ================================================================

int main (int argc, char* argv[])
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s K lambda abszed\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	double lambda = atof(argv[2]);
	double abszed = atof(argv[3]);
	double beta = 2.0*K;

// #define ZEE
#ifdef ZEE
	double complex z = abszed * (cos(M_PI*lambda) + I*sin(M_PI*lambda));
	print_vee(beta, z);
#endif

	int n =  lambda;
	print_d(n, beta);
}
