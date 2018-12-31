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

// This performs the sum to get the Renyi-Gelfond-Parry invariant
// measure, but rotated by exp(i\pi\lambda)
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

// Same as above, but for general z.
double complex rhoz(double x, double beta, complex double z)
{
	x /= 0.5*beta;
	if (1.0 < x) return 0.0;

	double ob = 1.0/beta;
	double obn = 1.0;
	complex double zn = 1.0;

	// accumulated sum
	complex double rho = 0.0;

	double tn = 1.0;

	while (1.0e-18 < obn)
	{
		if (x < tn) {
			rho += obn * zn;
		}

		// compute 1/beta^N
		obn *= ob;

		// compute xform^N(1)
		tn *= beta;
		if (1.0 < tn) tn -= 1.0;

		// compute z^n;
		zn *= z;
	}

	return rho;
}

// xfer operator
void cxfer(double complex* dest, double complex* src, int npts, double beta)
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

void do_unitary(double beta, double lambda)
{

#define NPTS 3801
	double complex rho[NPTS];
	double sum = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		rho[i] = rhogu(x, beta, lambda);
#if MODULUS_PHASE
		double re = creal(rho[i]);
		double im = cimag(rho[i]);
		double rms = sqrt(re*re + im*im);
		double phase = atan2(im, re);
#endif
		sum += creal(rho[i]);
	}

#ifdef RENORMLIZE
	double norm = sum / (double) NPTS;
	for (int i=0; i< NPTS; i++)
	{
		rho[i] /= norm;
	}
#endif

	double complex trho[NPTS];
	cxfer(trho, rho, NPTS, beta);

	double co = cos(M_PI*lambda);
	double si = -sin(M_PI*lambda);

	// This graphs the sine and cosine density
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		double re = creal(rho[i]);
		double im = cimag(rho[i]);
		double tre = creal(trho[i]);
		double tim = cimag(trho[i]);
		double gre = re*co - im*si;
		double gim = im*co + re*si;
		printf("%d	%g	%g	%g	%g	%g	%g	%g	%g\n", i, x,
			re, im,
			tre, tim,
			tre-gre, tim-gim, tre-gre-(tim-gim)
			);
	}
}

// Alsmot same as above, but for fixed x.
void do_const(double beta, double lambda)
{

#define NPTS 3801
	double complex rho[NPTS];
	double sum = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		rho[i] = rhogu(x, beta, lambda);
		sum += creal(rho[i]);
	}

	double complex trho[NPTS];
	cxfer(trho, rho, NPTS, beta);

	double co = cos(M_PI*lambda);
	double si = -sin(M_PI*lambda);

	int i = 0.2*NPTS;
	{
		double re = creal(rho[i]);
		double im = cimag(rho[i]);
		double tre = creal(trho[i]);
		double tim = cimag(trho[i]);
		double gre = re*co - im*si;
		double gim = im*co + re*si;
		printf("%d	%g	%g	%g	%g	%g	%g	%g	%g\n", i, lambda,
			re, im,
			tre, tim,
			tre-gre, tim-gim,
         sqrt((tre-gre)*(tre-gre)+(tim-gim)*(tim-gim))
			);
	}
}

void do_zee(double beta, double complex zee)
{

// #define NPTS 3801
	double complex rho[NPTS];
	double sum = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		rho[i] = rhoz(x, beta, zee);
		sum += creal(rho[i]);
	}

	double complex trho[NPTS];
	cxfer(trho, rho, NPTS, beta);

	// This graphs the sine and cosine density
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double)i + 0.5) / ((double) NPTS);
		double re = creal(rho[i]);
		double im = cimag(rho[i]);
		double tre = creal(trho[i]);
		double tim = cimag(trho[i]);
		double complex gee = rho[i] / zee;
		double gre = creal(gee);
		double gim = cimag(gee);
		printf("%d	%g	%g	%g	%g	%g	%g	%g	%g\n", i, x,
			re, im,
			tre, tim,
			tre-gre, tim-gim, tre-gre-(tim-gim)
			);
	}
}

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

#ifdef UNITARY_ONLY
	do_unitary(beta, lambda);
#endif // NITARY_ONLY

#ifdef ZEE
	double complex z = abszed * (cos(M_PI*lambda) + I*sin(M_PI*lambda));
	do_zee(beta, z);
#endif

#define EXPLORE_CONST
#ifdef EXPLORE_CONST
	// use with deangle.gplot
#define NLAM 501
	printf("#\n# Angular dependence for beta=%g\n#\n", beta);
	for (int i=0; i< NLAM; i++)
	{
		double lam = ((double)i + 0.5) / ((double) NLAM);
		do_const(beta, lam);
	}
#endif
}
