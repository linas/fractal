
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
#include <gsl/gsl_sf_zeta.h>

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
		double ct = cos(theta);
		double st = -sin(theta);
		double term = pow(2.0*M_PI*((double) n), -re_s);
		double rearg = cos (2.0*M_PI*((double) n)*x);
		double imarg = sin (2.0*M_PI*((double) n)*x);
		reacc += term * (rearg * ct - imarg*st);
		imacc += term * (rearg * st + imarg*ct);
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

#define DO_RENORM
#ifdef DO_RENORM
	double r = exp (lnr.val);
	reacc *= 2.0 *r;
	imacc *= 2.0 *r;
#endif

	double gc = cos (garg.val);
	double gs = sin (garg.val);

	double tmp = reacc *gc - imacc * gs;
	imacc = reacc *gs + imacc * gc;
	reacc = tmp;

	*rea = reacc;
	*ima = imacc;
}


/* with huge eignvalues, must take s less than zero only, though */

double beta_even_grow (double x, double sre)
{
	/* Hurwitz Zeta Function
	 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
	 *
	 * s > 1.0, q > 0.0
	 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
	 */

	double even = gsl_sf_hzeta (1.0-sre, x) + gsl_sf_hzeta (1.0-sre, 1.0-x);
	even *= sre;
	even /= cos (0.5*M_PI*sre);
	return even;
}

double beta_odd_grow (double x, double sre)
{
	double odd = gsl_sf_hzeta (1.0-sre, x) - gsl_sf_hzeta (1.0-sre, 1.0-x);
	odd *= sre;
	odd /= sin (0.5*M_PI*sre);
	return odd;
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

	sre = -s;
	sim = 0.0;

	// double lambda = pow (0.5, s);
	double lambda = pow (0.5, sre);
	double re_lambda = lambda * cos (sim*log(2.0));
	double im_lambda = - lambda * sin (sim*log(2.0));
	lambda = sqrt (re_lambda*re_lambda + im_lambda*im_lambda);

	printf ("#\n# ess=%g +i %g  eigenvalue lambda=%g +i %g \n#\n", sre, sim, re_lambda, im_lambda);
	printf ("#\n#   |lambda|=%g  arg(lambda)=%g \n#\n", lambda, sim*log(2.0));
	
	int imax = 623;
	for (i=1; i<imax; i++) 
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

#ifdef GRAPH_EIGENSTATES
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
#endif

// #define CHECK_THAT_ZETA_GROW_IS_EIGENSTATE 1
#if CHECK_THAT_ZETA_GROW_IS_EIGENSTATE 
		// verify its an eigenstate
		// double y = beta_even_grow (x,sre);
		// double z = 0.5*(beta_even_grow (0.5*x, sre) + beta_even_grow (0.5+0.5*x, sre));
		double y = beta_odd_grow (x,sre);
		double z = 0.5*(beta_odd_grow (0.5*x, sre) + beta_odd_grow (0.5+0.5*x, sre));
		z /= y;
		z /= lambda;
		z -= 1.0;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
#endif

#if SING_FREE_TYPE_ONE
		double e = beta_even_grow (x,sre);
		double o = beta_odd_grow (x,sre);

		e *= cos (0.5*M_PI*sre);
		o *= sin (0.5*M_PI*sre);
		double poz = pow (x, sre-1.0);
		double poo = pow (1.0-x, sre-1.0);
		e /= poz + poo;
		o /= poz - poo;
		e /= sre;
		o /= sre;

		e -= 1.0;
		o -= 1.0;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, e, o);
#endif

		double e = beta_even_grow (x,sre);
		double o = beta_odd_grow (x,sre);

		e *= cos (0.5*M_PI*sre);
		o *= sin (0.5*M_PI*sre);
		e /= sre;
		o /= sre;

		double poz = pow (x, sre-1.0);
		double poo = pow (1.0-x, sre-1.0);
		e -= poz + poo;
		o -= poz - poo;

		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, e, o);

		fflush (stdout);
	}
}
