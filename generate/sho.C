/*
 * sho.c
 *
 * High-precison Simple Harmonic Oscillator, using the 
 * Gnu Multiple-precision library.
 *
 * Actually, high-precision eigenfunctions for arbitrary eigenvalue
 * 
 * Linas Vepstas November 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "../misc/bignum/zmp.h"

/* ======================================================================= */
/* rough count of number of digits in a number */

static inline unsigned int num_digits (mpz_t num, mpz_t tmpa, mpz_t tmpb)
{
	unsigned int n=0;
	
	mpz_set (tmpb, num);
	while (1)
	{
		mpz_fdiv_q_ui (tmpa, tmpb, 100);
		mpz_set (tmpb, tmpa);
		if (0 == mpz_sgn  (tmpa)) break;
		n += 2;
	}
	return n;
}

/* ======================================================================= */
/* compute entire_sub_s for complex-valued s
 */

void psi_one (mpf_t re_psi, mpf_t im_psi, 
              double re_y, double im_y, unsigned int prec, int nterms)
{
	double reysq, imysq, tmp;

	double relam = 1.5;
	double imlam = 0.0;
	
	/* first, take the square of y */
	reysq = re_y*re_y - im_y*im_y;
	imysq = 2.0*re_y*im_y;
	
	mpf_set_ui (re_psi, 0);
	mpf_set_ui (im_psi, 0);

	cpx_t em, a, b;
	cpx_init (&em);
	cpx_init (&a);
	cpx_init (&b);

	cpx_confluent (&em, 
						 
	cpx_clear (&em);
	cpx_clear (&a);
	cpx_clear (&b);

}

/* ============================================================================= */

static mpf_t re_a, im_a;
static int prec;
static int nterms;

static void psi_init (void)
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
	
static double psi (double re_q, double im_q)
{
	double deno = re_q - 1.0;
	deno = deno*deno + im_q*im_q;
	deno = 1.0/deno;
	double re_s = 2.0*im_q* deno;
	double im_s = (re_q*re_q + im_q*im_q - 1.0) * deno;;

	// re_s = re_q;
	// im_s = im_q;

	psi_one (re_a, im_a, re_s, im_s, prec, nterms);

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
	psi_init();

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      // if (i%10==0) printf(" start row %d\n", i);
      printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = psi (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
