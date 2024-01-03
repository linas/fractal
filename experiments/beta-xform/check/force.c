/*
 * force.c
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

/* Some magic numbers!
 *
 * #define NHIST 27717  with "SAW" intialization and with beta=1.61803
 * converges in a stable way to the first decaying golden eigenfunc.
 * Other intializations are NOT stable and do NOT converge!
 * Other numbers of bins do NOT result in a stable convergence!
 * Even the explicit initialization of the eigenfunc is NOT  stable!
 * Other beta values are NOT stable, including beta=1.618034
 *
 * What is causing (in)stability?
 *
 * Other magic numbers:
 * Using the SAW or GOLD_ONE initialization:
 * #define NHIST 28658 (fibo 22) + 1 with beta=1.618034 has ringing pattern
 * #define NHIST any fib +1 seems to work, with ringing.
 */

// #define NHIST 360361
// #define NHIST 160361
// #define NHIST 27721
// #define NHIST 27719
// #define NHIST 27717
// #define NHIST 28657 // (fibo 22)
// #define NHIST 28658
// #define NHIST 196418 // (fib 26)
#define NHIST 196419
double histo[NHIST];
double histn[NHIST];

int enx (double x)
{
	int n = floor (NHIST * x);
	if (0 > n) n=0;
	if (NHIST <= n) n=NHIST-1;
	return n;
}

double dense(double x)
{
	return histo[enx(x)];
}

double ell(double beta, double y)
{
	if (0.5 * beta < y) return 0.0;

	double xlo = y / beta;
	double xhi = xlo + 0.5;

	double dlo = dense(xlo);
	double dhi = dense(xhi);

	double ellie = dlo + dhi;
	ellie /= beta;
	return ellie;
}

void setup(double beta)
{
#define SAW
#ifdef SAW
	for (int i=0; i<NHIST; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
	}
#endif

// #define STEP
#ifdef STEP
	int midp = enx(0.25*beta);
	for (int i=0; i<midp; i++) histn[i] = 1.0;
	for (int i=midp; i<2*midp; i++) histn[i] = -1.0;
	for (int i=2*midp; i<NHIST; i++) histn[i] = 0.0;
#endif

// #define GOLD_ONE
#ifdef GOLD_ONE
	int half = NHIST/2;
	for (int i=0; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*x-0.5;
	}
	int endp = enx(0.5*beta);
	for (int i=half; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
	}
	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;
#endif
}

double normalize(double beta)
{
	double delta =1.0 / ((double) NHIST);

	// Remove constant
	double sum = 0.0;
	for (int i=0; i<NHIST; i++) sum += histn[i];
	for (int i=0; i<NHIST; i++) histn[i] -= sum * delta;

	// Normalize
	double norm = 0.0;
	for (int i=0; i<NHIST; i++) norm += fabs(histn[i]);
	norm *= delta;
	for (int i=0; i<NHIST; i++) histo[i] = histn[i] / norm;

	return norm;
}

double step(double beta)
{
	double delta =1.0 / ((double) NHIST);

	// Run step
	for (int i=0; i<NHIST; i++)
	{
		double x = (((double) i) + 0.5) * delta;
		histn[i] = ell(beta, x);
	}

	return normalize(beta);
}

#define NCAP 10
double capt[NCAP][NHIST];

void capture(int n)
{
	for (int i=0; i<NHIST; i++) capt[n][i] = histo[i];
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta nsteps\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int nsteps = atoi(argv[2]);

	setup(beta);
	normalize(beta);

#define SHOW_NORMS
#ifdef SHOW_NORMS
	printf("#\n# beta=%g NHIST=%d\n#\n", beta, NHIST);
	for (int i=0; i< nsteps; i++)
	{
		double lam = step(beta);
		printf("%d	%g\n", i, lam);
	}
#endif

// #define SHOW_DENS
#ifdef SHOW_DENS
	printf("#\n# beta=%g\n#\n", beta);

	for (int i=0; i< nsteps; i++)
		step(beta);

	for (int i=0; i< NCAP; i++)
	{
		step(beta);
		capture(i);
	}

	for (int j=0; j< NHIST; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NHIST);
		printf("%d	%g", j, x);

		for (int i=0; i< NCAP; i++)
		{
			printf("	%g", capt[i][j]);
		}
		printf("\n");
	}
#endif
}
