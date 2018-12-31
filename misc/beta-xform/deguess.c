/*
 * deguess.c
 * Wild-ass guess, using a spining version of Parry-Gel'fond expression.
 *
 * Dec 2018
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double complex rhogu(double x, double beta, double lambda)
{
	x /= 0.5*beta;
	if (1.0 < x) return 0.0;

	double co = cos(M_PI*lambda);
	double si = sin(M_PI*lambda);
	double cn = 1.0;  // cos(N*pi*lambda)
	double sn = 0.0;
	double ob = 1.0/beta;
	double obn = 1.0;

	// accumulated sum
	double rho_r = 0.0;
	double rho_i = 0.0;

	double tn = 1.0;

	while (1.0e-18 < obn)
	{
		if (x < tn) {
			rho_r += obn * cn;
			rho_i += obn * sn;
		}

		// compute 1/beta^N
		obn *= ob;

		// compute xform^N(1)
		tn *= beta;
		if (1.0 < tn) tn -= 1.0;

		// compute cosine(N*lambda) sin (N*lambda)
		double tmp = co*cn - si*sn;
		sn = co*sn + si*cn;
		cn = tmp;
	}

	return rho_r + I* rho_i;
}

// Real-part-only xfer operator
void xfer(double* dest, double* src, int npts, double beta)
{
	for (int i=0; i<npts; i++)
	{
		double y = ((double)i + 0.5) / ((double) npts);
		if (0.5*beta < y)
		{
			dest[i] = 0.0;
		}
		else
		{
			double x = npts * y / beta;
			int bin = floor(x);
			dest[i] = (src[bin] + src[bin+npts/2]) / beta;
		}
	}
}

int main (int argc, char* argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s K lambda\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	double lambda = atof(argv[2]);
	double beta = 2.0*K;

#define NPTS 1801
	double rerho[NPTS];
	double imrho[NPTS];
	double sum = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		double complex rho = rhogu(x, beta, lambda);
		double re = creal(rho);
		double im = cimag(rho);
		rerho[i] = re;
		imrho[i] = im;
#if INSANE
		double rms = sqrt(re*re + im*im);
		if (1.0e-14 < rms)
		{
			rerho[i] = re / rms;
			imrho[i] = im / rms;
		}
#endif
		sum += rerho[i];
	}
	double norm = sum / (double) NPTS;
	for (int i=0; i< NPTS; i++)
	{
		rerho[i] /= norm;
		imrho[i] /= norm;
	}

	double trerho[NPTS];
	xfer(trerho, rerho, NPTS, beta);

	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		printf("%d	%g	%g	%g	%g\n", i, x, rerho[i], imrho[i], trerho[i]);
	}
}
