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

	while (1.0e-16 < obn)
	{
		if (x < tn) {
			rho_r += obn * cn * tn;
			rho_i += obn * sn * tn;
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


int main (int argc, char* argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s K lambda\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	double lambda = atof(argv[2]);
	double beta = 2.0*K;

#define NPTS 801
	double rerho[NPTS];
	double imrho[NPTS];
	double ms = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		double complex rho = rhogu(x, beta, lambda);
		rerho[i] = creal(rho);
		imrho[i] = cimag(rho);
		ms += rerho[i]*rerho[i] + imrho[i]*imrho[i];
	}
	double rms = sqrt(ms / (double) NPTS);
	for (int i=0; i< NPTS; i++)
	{
		rerho[i] / = rms;
		imrho[i] / = rms;
	}
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		printf("%d	%g	%g	%g\n", i, x, rerho[i], imrho[i]);
	}

}
