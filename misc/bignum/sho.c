
/*
 * Simple Harmonic Oscilltor
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-hyper.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* ======================================================== */
/* The two eigenfunctions of the simple harmonic oscillator */

void psi_1 (cpx_t psi, cpx_t lambda, cpx_t y, int prec)
{
	cpx_t a, b, z;
	cpx_init (a);
	cpx_init (b);
	cpx_init (z);

	/* b=1/2 */
	cpx_set_ui (b, 1, 0);
	mpf_div_ui (b[0].re, b[0].re, 2);

	/* a= 1/4-lambda/2 */
	cpx_set (a,b);
	cpx_sub (a, a, lambda);
	cpx_div_ui (a, a, 2);

	/* z = y^2 */
	cpx_mul (z, y, y);

	cpx_confluent (psi, a, b, z, prec);

	/* psi_1 = exp (-y^2/2) * M(a,b,z) */
	cpx_div_ui (z, z, 2);
	cpx_neg (z, z);
	cpx_exp (z, z, prec);
	cpx_mul (psi, psi, z);
	
	cpx_clear (a);
	cpx_clear (b);
	cpx_clear (z);
}

void psi_2 (cpx_t psiret, cpx_t lambda, cpx_t y, int prec)
{
	cpx_t a, b, z, psi;
	cpx_init (a);
	cpx_init (b);
	cpx_init (z);
	cpx_init (psi);

	/* b=3/2 */
	cpx_set_ui (b, 3, 0);
	cpx_div_ui (b, b, 2);

	/* a= (3/2-lambda)/2 */
	cpx_set (a,b);
	cpx_sub (a, a, lambda);
	cpx_div_ui (a, a, 2);

	/* z = y^2 */
	cpx_mul (z, y, y);

	cpx_confluent (psi, a, b, z, prec);

	/* psi_2 = y * exp (-y^2/2) * M(a,b,z) */
	cpx_div_ui (z, z, 2);
	cpx_neg (z, z);
	cpx_exp (z, z, prec);
	cpx_mul (psi, psi, z);
	cpx_mul (psi, psi, y);
	
	/* copy now, thus avoiding possible clobber of args */
	cpx_set (psiret, psi);

	cpx_clear (a);
	cpx_clear (b);
	cpx_clear (z);
	cpx_clear (psi);
}

void eta_1 (cpx_t eta, cpx_t lam, cpx_t y, int prec)
{
	cpx_t pha, lambda, l2;
	cpx_init (pha);
	cpx_init (lambda);
	cpx_init (l2);

	/* Make copy of argument NOW! */
	cpx_set (lambda, lam);
	cpx_div_ui (l2, lambda, 2);
	
	psi_1 (eta, lambda, y, prec);

	mpf_t tmp;
	mpf_init (tmp);

	/* phase term, exp (i pi lambda/2) */
	fp_pi (tmp, prec);
	cpx_mul_mpf (pha, l2, tmp);
	cpx_times_i (pha, pha);
	cpx_exp (pha, pha, prec);
	cpx_mul (eta, eta, pha);

	/* times sqrt(2pi) */
	fp_sqrt_two_pi (tmp, prec);
	cpx_mul_mpf(eta, eta, tmp);
	
	/* power of two term */
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 2);
	cpx_pow_mpf (pha, tmp, l2, prec);
	cpx_mul (eta, eta, pha);
	
	/* Gamma (3/4+lambda/2) */
	mpf_set_ui (tmp, 3);
	mpf_div_ui (tmp, tmp, 4);
	cpx_add_mpf (pha, l2, tmp);
	cpx_gamma (pha, pha, prec);
	cpx_div (eta, eta, pha);

	cpx_clear (pha);
	cpx_clear (lambda);
	cpx_clear (l2);
	mpf_clear (tmp);
}

void eta_2 (cpx_t eta, cpx_t lam, cpx_t y, int prec)
{
	cpx_t pha, lambda, l2;
	cpx_init (pha);
	cpx_init (lambda);
	cpx_init (l2);

	/* Make copy of arg NOW! */
	cpx_set (lambda, lam);
	cpx_set (l2, lam);
	cpx_div_ui (l2, l2, 2);

	psi_2 (eta, lambda, y, prec);

	mpf_t tmp;
	mpf_init (tmp);

	/* phase term, exp (i pi (lambda-1)/2) */
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 2);
	cpx_sub_mpf (pha, l2, tmp);
	fp_pi (tmp, prec);
	cpx_mul_mpf (pha, pha, tmp);
	cpx_times_i (pha, pha);
	cpx_exp (pha, pha, prec);
	cpx_mul (eta, eta, pha);

	/* times sqrt(2pi) */
	fp_sqrt_two_pi (tmp, prec);
	cpx_mul_mpf(eta, eta, tmp);
	
	/* power of two term */
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 2);
	cpx_pow_mpf (pha, tmp, l2, prec);
	cpx_mul (eta, eta, pha);
	cpx_mul_ui (eta, eta, 2);
	
	/* Gamma (1/4+lambda/2) */
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 4);
	cpx_add_mpf (pha, l2, tmp);
	cpx_gamma (pha, pha, prec);
	cpx_div (eta, eta, pha);

	cpx_clear (pha);
	cpx_clear (lambda);
	cpx_clear (l2);
	mpf_clear (tmp);
}

/* ======================================================== */

void validate_ratio (cpx_t lambda, cpx_t y, int prec)
{
	cpx_t e1, e2, lam, rat, tang;
	cpx_init (e1);
	cpx_init (e2);
	cpx_init (lam);
	cpx_init (rat);
	cpx_init (tang);

	mpf_t pi;
	mpf_init (pi);

	int n;
	for (n=4; n<14; n+=3)
	{
		printf ("n=%d  \n", n);

		/* lambda-n */
		cpx_add_ui (lam, lambda, n, 0);
		
// #define VALIDATE_PSI_ONE
#ifdef VALIDATE_PSI_ONE
		/* The following validates just great, 
		 * for both large positive and large negative lambda */
		psi_1 (e1, lam, y, prec);
		cpx_prt ("psi_1=", e1);
		printf ("\n");

		double flam = mpf_get_d (lam[0].re);
		double fy = mpf_get_d (y[0].re);
		double s1 = sqrt (2.0*flam-1);
		double cs1 = cos (fy*s1);

		double sn1 = - 0.5*fy*sin(fy*s1)*(1.0-fy*fy/3.0)/s1;
		// cs1 += sn1;
		// double cs1 = 0.5* exp (fy*s1);
		double fpsi1 = mpf_get_d (e1[0].re);
		printf ("duude cs1=%g sn1=%g sum=%g cmp=%g\n", cs1, sn1,cs1+sn1, (cs1+sn1)/fpsi1);
#endif

// #define VALIDATE_PSI_TWO
#ifdef VALIDATE_PSI_TWO
		/* The following validates just great, for large negative
		 * lambda, and seems numerically unstable for pos lambda */
		psi_2 (e2, lam, y, prec);
		cpx_prt ("e2=", e2);
		printf ("\n");

		double flam = mpf_get_d (lam[0].re);
		double fy = mpf_get_d (y[0].re);
		double s2 = sqrt (3.0-2.0*flam);
		// double cs2 = 0.5 * cos (fy*s2);
		double cs2 = 0.5 * exp (fy*s2) / s2;
		double fe2 = mpf_get_d (e2[0].re);
		printf ("duude cs1=%g cmp=%g\n", cs2, cs2/fe2);
#endif

#if VALIDATE_PSI_RATIO
		/* The following validates just great, 
		 * for large negative lambda. */
		psi_1 (e1, lam, y, prec);
		psi_2 (e2, lam, y, prec);
		cpx_div (rat, e2, e1);
		cpx_prt ("rat=", rat);
		printf ("\n");

		double flam = mpf_get_d (lam[0].re);
		double fy = mpf_get_d (y[0].re);
		double s1 = sqrt (1.0-2.0*flam);
		double s2 = sqrt (3.0-2.0*flam);
		double frat = exp (fy*(s2-s1)) / s2;
		frat = 1.0 / s2;
		double fr = mpf_get_d (rat[0].re);

		printf ("duude frat=%g rr=%g\n", frat, frat/fr);
#endif

#if VALIDATE_ETA_RATIO
		eta_1 (e1, lam, y, prec);
		eta_2 (e2, lam, y, prec);
		cpx_div (rat, e2, e1);
		cpx_prt ("rat=",rat);
		printf ("\n");

		// double fr = mpf_get_d (rat[0].im);
		// fr /= tan (M_PI * 0.5*(mpf_get_d(lam[0].re)+0.5));
		// printf ("duude fr=%g\n", fr);
		
		/* compute lambda/2 + 1/4 */
		cpx_set_ui (tang, 1,0);
		cpx_div_ui (tang, tang, 2);
		cpx_add (tang, tang, lam);
		cpx_div_ui (tang, tang, 2);
		
		/* compute tan(lambda/2 + 1/4) */
		fp_pi (pi, prec);
		cpx_mul_mpf (tang, tang, pi);
		cpx_tangent (tang, tang, prec);
		cpx_times_i (tang, tang);
		
		cpx_prt ("tan=", tang);
		printf ("\n");
		
		/* ratio */
		cpx_div (rat, rat, tang);
		cpx_neg (rat, rat);
		cpx_sub_ui (rat, rat, 1, 0);


#endif

		printf ("\n");
	}

	mpf_clear (pi);
	cpx_clear (e1);
	cpx_clear (e2);	
	cpx_clear (lam);
	cpx_clear (rat);
	cpx_clear (tang);
}

/* ======================================================== */
/* Validate the asymptotic expansion for the gamma */
#if 0
void wtf (void)
{
	double lam;

	/* the plus-infty direction */
	for (lam=2.321; lam<30; lam += 0.3)
	{
		double r = tgamma (0.5*lam-0.25);
		r /= tgamma (0.5*lam+0.25);
		r *= lam-0.5;
		r /= sqrt (2*lam);

		printf ("its lam=%g\t rat=%g\n", lam, r);
	}
	
	/* The minus-infty direction, which is what we want */
	for (lam=-2.321; lam>-30; lam -= 0.3)
	{
		double r = tgamma (0.5*lam-0.25);
		r /= tgamma (0.5*lam+0.25);
		r *= lam-0.5;
		r /= sqrt (-2*lam);

		double a = tan (M_PI*(0.5*lam+0.25));

		printf ("its lam=%g\t rat=%g\t asymp=%g\t rr=%g\n", lam, r, a, r/a);
	}
}
#endif

/* ======================================================== */

void get_ratio (cpx_t tang, cpx_t lambda, int prec)
{
	mpf_t pi;
	mpf_init (pi);
	fp_pi (pi, prec);

	cpx_t lam;
	cpx_init (lam);
	cpx_set (lam, lambda);

	/* compute lambda/2 + 1/4 */
	cpx_set_ui (tang, 1,0);
	cpx_div_ui (tang, tang, 2);
	cpx_add (tang, tang, lam);
	cpx_div_ui (tang, tang, 2);
		
	/* compute tan(lambda/2 + 1/4) */
	cpx_mul_mpf (tang, tang, pi);
	cpx_tangent (tang, tang, prec);
	cpx_times_i (tang, tang);

	cpx_recip (tang, tang);

	cpx_clear (lam);
	mpf_clear (pi);
}

void coherent (cpx_t coho, cpx_t lambda, cpx_t que, cpx_t y, int prec)
{
	cpx_t e1, e2, rat, qn, lam;
	cpx_init (e1);
	cpx_init (e2);
	cpx_init (rat);
	cpx_init (qn);
	cpx_init (lam);

	cpx_set (qn, que);


	eta_1 (e1, lambda, y, prec);
	eta_2 (e2, lambda, y, prec);
	get_ratio (rat, lambda, prec);
	cpx_mul (e2,e2,rat);
	cpx_sub (coho, e1, e2);
	
	int n;
	for (n=1; n<140; n++)
	{
	printf ("n=%d  \n", n);
		/* lambda+n */
		cpx_add_ui (lam, lambda, n, 0);
		eta_1 (e1, lam, y, prec);
		eta_2 (e2, lam, y, prec);

		get_ratio (rat, lam, prec);
		cpx_mul (e2,e2,rat);
		cpx_sub (e1, e1, e2);

		cpx_mul (e1, e1, qn);
		cpx_add (coho, coho, e1);
		
		/* lambda-n */
		cpx_sub_ui (lam, lambda, n, 0);
		eta_1 (e1, lam, y, prec);
		eta_2 (e2, lam, y, prec);

	cpx_prt ("neg e1 =", e1);
	printf ("\n");
		get_ratio (rat, lam, prec);
		cpx_mul (e2,e2,rat);
		cpx_sub (e1, e1, e2);

	cpx_prt ("neg thing=", e1);
	printf ("\n");
		cpx_div (e1, e1, qn);
		cpx_add (coho, coho, e1);
		
	cpx_prt ("coho=", coho);
	printf ("\n");
	printf ("\n");
		
		cpx_mul (qn, qn, que);
	}

	cpx_clear (e1);
	cpx_clear (e2);	
	cpx_clear (rat);	
	cpx_clear (qn);
	cpx_clear (lam);
}

/* ======================================================== */

void prt_graph (double lam, int prec)
{
	cpx_t ps1, ps2, lambda, z;
	cpx_init (ps1);
	cpx_init (ps2);
	cpx_init (lambda);
	cpx_init (z);

	cpx_set_ui (lambda, 1, 0);
	mpf_set_d (lambda[0].re, lam);
	mpf_set_d (lambda[0].im, 0.0);

	printf ("# \n# Graph of simple harmonic oscillator\n# \n");
	cpx_prt ("# lambda = ", lambda);
	printf ("\n# \n");
	fflush (stdout);

	cpx_t ex;
	cpx_init (ex);
	cpx_set_ui (ex, 0,0);

	double x;
	for (x=-6.0; x<6.01; x+=0.1)
	{
		mpf_set_d (ex[0].re, x);
		mpf_set_d (ex[0].im, 1.0);

#if 0
		psi_1 (ps1, lambda, ex, prec);
		psi_2 (ps2, lambda, ex, prec);
#endif

		eta_1 (ps1, lambda, ex, prec);
		eta_2 (ps2, lambda, ex, prec);

		printf ("%g", x);
		// fp_prt ("\t", ps1[0].re);
		// fp_prt ("\t", ps1[0].im);
		// fp_prt ("\t", ps2[0].re);
		// fp_prt ("\t", ps2[0].im);

		/* phase */
		fp_arctan2 (ps1[0].re, ps1[0].im, ps1[0].re, 30);
		fp_arctan2 (ps2[0].re, ps2[0].im, ps2[0].re, 30);

		mpf_sub (ps1[0].re, ps1[0].re, ps2[0].re);

		fp_prt ("\t", ps1[0].re);
		fp_prt ("\t", ps2[0].re);
		printf ("\n");
		fflush (stdout);
	}

}

void do_coho (double lam, int prec)
{
	cpx_t coho, lambda, q, z;
	cpx_init (coho);
	cpx_init (lambda);
	cpx_init (q);
	cpx_init (z);

	cpx_set_d (lambda, lam, 0.2);

	cpx_set_d (q, 4.1, 0.0);
	cpx_set_d (z, 2.5671, 0.2);
	
	// coherent (coho, lambda, q, z, prec);
	validate_ratio (lambda, z, prec);
}

/* ======================================================== */

int
main (int argc, char *argv[])
{
	int prec = 50;
	
	if (2> argc)
	{
		fprintf (stderr, "Usage: %s lambda\n", argv[0]);
		exit (1);
	} 

	double lam;
	lam = atof (argv[1]);

	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);
	
	// prt_graph (lam, prec);
	do_coho (lam, prec);
	// do_coho (lam, prec);
	
	return 0;
}
