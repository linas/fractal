/*
 * reigen.c
 *
 * Recursve eigenfunctions for the undershift
 * Dec 2017
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Compute eignefunction, recursively.
double reig(double x, double K, double lambda, int niter)
{
	if (K < x) return 0.0;
	if (niter < 0)
	{
		// Approximate by a constant.
		if (0.999999999 < lambda) return 1.0 / K;

		// Approximate by something that integrates to zero.
		if (x < 0.25*K) return 1.0/K;
		if (x > 0.75*K) return 1.0/K;
		return -1.0/K;

		// Approximate by something that integrates to zero.
		if (x < 0.5*K) return 1.0/K;
		return -1.0/K;
	}

	double tkay = 2.0*K;

	// Short-cut. The other branch has vanishing contribution,
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, lambda, niter-1) / (lambda * tkay);
	}
	double sum = reig(x/tkay, K, lambda, niter-1);
	sum += reig(0.5 + x/tkay, K, lambda, niter-1);
	sum /= lambda * tkay;
	return sum;
}

// Apply the xfer operator. A non-recursive version of above.
void apply_xfer(double *out, double *in, int sz, double K)
{
	double tkay = 2.0 * K;
	for (int i=0; i<sz; i++)
	{
		double x = ((double) i + 0.5) / ((double) sz);
		if (x < K)
		{
			double b1 = x / tkay;
			double b2 = b1 + 0.5;
			unsigned int isamp1 = floor(b1 * sz);
			unsigned int isamp2 = floor(b2 * sz);
			if (isamp1 >= sz) printf ("ohh noo %d %g\n", isamp1, b1*sz);
			assert (isamp1< sz); 
			if (isamp2 >= sz) printf ("ohh 2 noo %d %g\n", isamp2, b2*sz);
			assert (isamp2< sz); 
			out[i] = (in[isamp1] + in [isamp2]) / tkay;
		}
		else
			out[i] = 0.0;
	}
}

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K w\n", argv[0]);
		exit (1);
	}
	double K = atof(argv[1]);
	double lambda = atof(argv[2]);

	printf("#\n# K = %g lambda = %g\n#\n", K, lambda);

#define NPTS 1803

	double psi[NPTS];

// #define NRECU 40
#define NRECU 20

	// Compute an eigenfunction, recursively.
	double lin = 0.0;
	double squ = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double y = reig(x, K, lambda, NRECU);
		lin += y;
		squ += y*y;

		psi[i] = y;
	}
	lin /= NPTS;
	squ = sqrt(squ / NPTS);

	printf ("# linear=%g squ=%g\n#\n", lin, squ);

	// Renormalize to unit mean-square 
	for (int i=0; i<NPTS; i++)
	{
		psi[i] /= squ;
	}

	// Apply the xfer operator.
	double verify[NPTS];
	apply_xfer(verify, psi, NPTS, K);

	// See how the computed and expected functions differ
	double sumerr = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double y = psi[i];
		double z = verify[i] / lambda;
		double err = y-z;
		sumerr += err*err;
	}
	sumerr /= NPTS;
	printf("# eeign error: %g\n#\n", sumerr);

	// Dump to file.
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double y = psi[i];
		double z = verify[i];
		printf("%d	%g %g	%g\n", i, x, y, z);
	}
}
