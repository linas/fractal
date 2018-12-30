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

double complex rho(double x, double beta, double lambda)
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
	double K = atof(argv[1]);
	double x = atof(argv[2]);
	double beta = 2.0*K;

}
