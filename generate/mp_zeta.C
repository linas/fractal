/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Actually, high-precision a_s on the complex plane
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "../misc/bignum/zmp.h"

/* ======================================================================= */
/* compute entire_sub_s for complex-valued s
 */

void entire_sub_s (mpf_t re_b, mpf_t im_b, double re_s, double im_s, unsigned int prec, int nterms)
{
	int k;
	mpf_t rebin, imbin, term, racc, iacc, rzeta, izeta;

	mpf_init (term);
	mpf_init (racc);
	mpf_init (iacc);
	mpf_init (rzeta);
	mpf_init (izeta);
	mpf_init (rebin);
	mpf_init (imbin);
	
	mpf_set_ui (re_b, 0);
	mpf_set_ui (im_b, 0);

	int n = 650;  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	n = nterms;
	int downer  = 0;
	for (k=2; k<= n; k++)
	{
		/* Commpute the binomial */
		c_binomial_d (rebin, imbin, re_s, im_s, k);

// printf ("duude s= (%g %g) k=%d bin=(%g %g)\n", re_s, im_s, k, mpf_get_d(rebin), mpf_get_d(imbin));

		/* compute zeta (k) */
		fp_zeta (term, k, prec);
		mpf_sub_ui (term, term, 1);

		mpf_mul (rzeta, term, rebin);
		mpf_mul (izeta, term, imbin);

		if (k%2)
		{ 
			mpf_sub (racc, re_b, rzeta);
			mpf_sub (iacc, im_b, izeta);
		}
		else 
		{
			mpf_add (racc, re_b, rzeta);
			mpf_add (iacc, im_b, izeta);
		}
		
		mpf_set (re_b, racc);
		mpf_set (im_b, iacc);
#if 1
		double rt = mpf_get_d (rzeta);
		double it = mpf_get_d (izeta);
		double ra = mpf_get_d (re_b);
		double ia = mpf_get_d (im_b);
		if (rt*rt +it*it < 1.0e-15 * (ra*ra+ia*ia)) 
		{
			if (downer > 5) break;
			downer ++;
		}
#endif

	}

	mpf_clear (term);
	mpf_clear (racc);
	mpf_clear (iacc);
	mpf_clear (rzeta);
	mpf_clear (izeta);
	mpf_clear (rebin);
	mpf_clear (imbin);
}

/* ============================================================================= */

static mpf_t re_a, im_a;
static int prec;
static int nterms;

static void a_s_init (void)
{
	/* the decimal precison (number of decimal places) */
	prec = 300;
   nterms = 300;

	prec = 100;
   nterms = 100;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 30);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_init (re_a);
	mpf_init (im_a);
}
	
static double a_s (double re_q, double im_q)
{
	double deno = re_q - 1.0;
	deno = deno*deno + im_q*im_q;
	deno = 1.0/deno;
	double re_s = 2.0*im_q* deno;
	double im_s = (re_q*re_q + im_q*im_q - 1.0) * deno;;

	// re_s = re_q;
	// im_s = im_q;

	// a_sub_s (re_a, im_a, re_s, im_s, prec);
	// b_sub_s (re_a, im_a, re_s, im_s, prec, nterms);
	entire_sub_s (re_a, im_a, re_s, im_s, prec, nterms);

	double frea = mpf_get_d (re_a);
	double fima = mpf_get_d (im_a);

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler totient in a simple way 
 */

void 
MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

	prec = itermax;
	a_s_init();

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      // if (i%10==0) printf(" start row %d\n", i);
      printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = a_s (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
