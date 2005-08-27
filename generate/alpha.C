/*
 * alpha.C
 *
 * FUNCTION:
 * display the alpha(z) function from the reiman-pochhammer  paper.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>

#include "brat.h"

static void alpha_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	/*  1-z */
	double re_one_minus_z = 1.0 - re_q;
	double im_one_minus_z = - re_q;

	/* 1/(1-z) */
	tmp = re_one_minus_z*re_one_minus_z + im_one_minus_z*im_one_minus_z;
	tmp = 1.0/tmp;
	double re_one_over_one_minus_z = re_one_minus_z * tmp;
	double im_one_over_one_minus_z = -im_one_minus_z * tmp;

	/* 1/z */
	tmp = 1.0/(re_q*re_q + im_q * im_q);
	double re_one_over_z = re_q * tmp;
	double im_one_over_z = -im_q * tmp;
			  
	/* ln gamma (1/(1-z)) */
	gsl_sf_result mod_lng, arg_lng;
	gsl_sf_lngamma_complex_e (re_one_over_one_minus_z, im_one_over_one_minus_z, 
	                &mod_lng, &arg_lng);

	double re_lng = mod_lng.val * cos (arg_lng.val);
	double im_lng = mod_lng.val * sin (arg_lng.val);

	/* (1/z) * ln gamma (1/(1-z)) */
	rep += re_one_over_z * re_lng - im_one_over_z * im_lng;
	imp += re_one_over_z * im_lng + im_one_over_z * re_lng;

	/* log (1-z) */
	gsl_sf_result m_lnz, t_lnz;
	gsl_sf_complex_log_e (re_one_minus_z, im_one_minus_z, &m_lnz, &t_lnz);
	double re_ln_omz = m_lnz.val;
	double im_ln_omz = t_lnz.val;
	
	/* (1/z) * log (1-z) */
	double re_log_over_z = re_ln_omz * re_one_over_z - im_ln_omz * im_one_over_z;
	double im_log_over_z = re_ln_omz * im_one_over_z + im_ln_omz * re_one_over_z;

	/* add (1/z) * log (1-z) * (1/(1-z)) */
	rep += re_log_over_z * re_one_over_one_minus_z - im_log_over_z * im_one_over_one_minus_z;
	imp += re_log_over_z * im_one_over_one_minus_z + im_log_over_z * re_one_over_one_minus_z;
	
	/* add - 0.5 * (1/z) * log (1-z) */
	rep += -0.5 * re_log_over_z;
	imp += -0.5 * im_log_over_z;

	/* add 1/(1-z) */
	rep += re_one_over_one_minus_z;
	imp += im_one_over_one_minus_z;
	
	*prep = rep;
	*pimp = imp;
}

static double alpha_series (double re_q, double im_q)
{
	double rep, imp;
	alpha_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return imp;
	return (atan2 (imp, rep) + M_PI) / (2.0*M_PI);
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler alpha in a simple way 
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

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = alpha_series (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
