
/* 
 * zeta.C
 * 
 * Explore Hurwitz zeta eigenfunctions of Bernoulli map
 *
 * Linas November 2004
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>

/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
// int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


double beta (double x, double s)
{
	int n;
	double acc = 0.0;
	for (n=1; n<57550; n++)
	{
		double term = pow(2.0*M_PI*((double) n), -s);
		acc += term * cos (2.0*M_PI*((double) n)*x);
		// acc += term * sin (2.0*M_PI*((double) n)*x);
		if (1.0e-16>term) break;
	}
	acc *= 2.0 * tgamma (s+1.0);
	return acc;
}

void betas (double x, double re_s, double im_s, double *rea, double *ima)
{
	int n;
	double reacc = 0.0;
	double imacc = 0.0;
	for (n=1; n<57550; n++)
	{
		double theta = im_s * log (2.0*M_PI*((double) n));
		double term = pow(2.0*M_PI*((double) n), -re_s);
		double arg = cos (2.0*M_PI*((double) n)*x);
		reacc += term * arg * cos(theta);
		imacc += -term * arg * sin(theta);
		// acc += term * sin (2.0*M_PI*((double) n)*x);
		if (1.0e-16>term) break;
	}

	// Urgh. Now normalize properly
	// acc *= 2.0 * tgamma (s+1.0);

	gsl_sf_result lnr, garg;
	int err = gsl_sf_lngamma_complex_e (re_s+1.0, im_s, &lnr, &garg);
	if (err)
	{
		printf ("# uhh ohhhh ! s = %g + i %g\n", re_s, im_s); 
	}
	// double r = exp (lnr.val);
	// reacc *= 2.0 *r;
	// imacc *= 2.0 *r;

	double gc = cos (garg.val);
	double gs = sin (garg.val);

	double tmp = reacc *gc - imacc * gs;
	imacc = reacc *gs + imacc * gc;
	reacc = tmp;

	*rea = reacc;
	*ima = imacc;
}

main (int argc, char * argv[])
{
	int i;

	double s=3.345;

	if (2>argc)
	{
		printf ("Usage: %s  <s-value>\n", argv[0]);
		exit (1);
	}
	s = atof (argv[1]);

	double sre=2.345;
	double sim = s;

	// double lambda = pow (0.5, s);
	double lambda = pow (0.5, sre);
	double re_lambda = lambda * cos (sim*log(2.0));
	double im_lambda = - lambda * sin (sim*log(2.0));
	lambda = sqrt (re_lambda*re_lambda + im_lambda*im_lambda);

	printf ("#\n# ess=%g +i %g  eigenvalue lambda=%g +i %g \n#\n", sre, sim, re_lambda, im_lambda);
	printf ("#\n#   |lambda|=%g  arg(lambda)=%g \n#\n", lambda, sim*log(2.0));
	
	int imax = 153;
	for (i=0; i<=imax; i++) 
	{
		double x = i/((double) imax);
		// double y = beta (x,s);

// #define CHECK_THAT_ZETA_IS_EIGENSTATE 1
#if CHECK_THAT_ZETA_IS_EIGENSTATE 
		// verify its an eigenstate
		double y = beta (x,s);
		double z = 0.5*(beta (0.5*x, s) + beta (0.5+0.5*x, s));
		z /= lambda;
		z -= y;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
#endif

// #define CHECK_COMPLEX_EIGENSATES
#ifdef CHECK_COMPLEX_EIGENSATES
		double rey, imy;
		betas (x, sre, sim, &rey, &imy);

		double reb1, imb1, reb2, imb2;

		betas (0.5*x, sre, sim, &reb1, &imb1);
		betas (0.5+0.5*x, sre, sim, &reb2, &imb2);

		double zre = 0.5*(reb1 + reb2);
		double zim = 0.5*(imb1 + imb2);

		double reyl = rey*re_lambda - imy*im_lambda;
		double imyl = rey*im_lambda + imy*re_lambda;

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, rey, zre-reyl, zim-imyl);
#endif

#define NPTS 5
		double rey, imy;
		int k;

		printf ("%d	%8.6g	", i, x);
		double esim = sim;
		for (k=0; k<NPTS; k++)
		{
			betas (x, sre, esim, &rey, &imy);
			// if (1==k%2) { rey = -rey; imy = -imy; }
			printf ("%8.6g	%8.6g	", rey, imy);
			esim += 2.0*M_PI / log(2.0);
		}
		printf ("\n");
		fflush (stdout);
	}
}
