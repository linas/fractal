/*
 * recheck.c
 * Recheck old (negative) results.
 * They're still negative.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
	return sum;
}

#define NHIST 1000
double histo[NHIST];

double nhits = 0.0;

double phi = 1.618033988749895;
double phi2 = 0.8090169943749475;
double zig(double x)
{
	x -= floor(x);
	if (x < phi2) return x/phi2;
	return (x-phi2) / (1-phi2);
}
double zigzag(double x)
{
	double w = 0.4;
	double tn = 1.0;
	double wn = 1.0;
	double sum = 0;
	for (int i=0; i<10; i++)
	{
		sum += wn * zig(tn*x);
		tn *= 2;
		wn *= w;
	}
	return sum;
}

double base(double x, double beta, double strength)
{
	int n = floor(NHIST * x);
	histo[n] += strength;
	nhits += strength;

	// double y = x * 2.0 / beta;
	// double ska = 1.0 / invar(beta, x);
	// return 1.0;
	// return 1.0 / invar(beta, x);

	// return 0.5-x;
	return x-0.5;
	// return 1.0-x;
	// return (y-0.5) *ska;

	// return sin(2.0* M_PI * x);
	// return sin(2.0* M_PI * y) / invar(beta, x);
	// return sin(2.0* M_PI * y) * ska;

	// return sqrt(invar(beta, x)) * ska;

#if ZIGZAG
	double a = 0.2;
	double b = 1.4473 / 1.61803;
	if (x< 0.5)
		return b + a*zigzag(2.0*x);
	if (x < 0.5*phi)
		return b+a*zigzag((x-0.5) * phi * 2);
	return 0.0;
#endif

	// if (x < 0.8*0.8) return 0.0;
// printf("yo x=%g\n", x);
	// if (x < 0.3333) return 0.0;
	// return 1.0;
}

double dense(double beta, double x)
{
	double midpnt = 0.5*beta;
	double obn = 1.0;
	double sum = 0.0;
	for (int i=0; i<1000; i++)
	{
		if (x < midpnt)
		{
			sum += obn * base (x, beta, obn);
		}
		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		// printf("i=%d x=%g midpnt=%g obn=%g sum=%g\n", i, x, midpnt, obn, sum);

		if (obn < 1e-15) break;
	}
	return sum;
}

double ell(double beta, double y)
{
	if (0.5 * beta < y) return 0.0;

	double xlo = y / beta;
	double xhi = xlo + 0.5;

	double dlo = dense(beta, xlo);
	double dhi = dense(beta, xhi);

	double ellie = dlo + dhi;
	ellie /= beta;
	return ellie;
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta \n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);

	for (int i=0; i<NHIST; i++) histo[i] = 0.0;

	int npts = 1000;
	printf("#\n# beta=%g\n#\n", beta);
	for (int i=0; i< npts; i++)
	{
		double x = (((double) i) + 0.5) / ((double) npts);
		dense(beta, x);
#if FAIL
		double f = dense(beta, x);
		double elf = ell(beta, x);
		double elg = ell(beta, beta*x);
		printf("%d	%g	%g	%g	%g\n", i, x, f, elf, elg);
#endif
	}

	double norm = nhits / NHIST;
	double save_histo[NHIST];
	for (int i=0; i< NHIST; i++) save_histo[i] = histo[i];

	printf("# norm = %g\n", norm);
	for (int i=0; i< NHIST; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		double h = save_histo[i] / norm;
		double f = dense(beta, x) / norm;
		double mu = invar(beta, x);
		double elf = ell(beta, x) / norm;
		printf("%d	%g	%g	%g	%g	%g\n", i, x, h, f, mu, elf);
	}
}
