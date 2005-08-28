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

static void const_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	int i;
	double tmp;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	double a[8];
	a[0] = -0.077215664901532851;
	a[1] = -0.0047486314774196408;
	a[2] = 0.00036610089349548089;
	a[3] = 0.00037600730566372326;
	a[4] = 0.00014301182486231440;
	a[5] = 3.3997818936021684e-5;
	a[6] = -4.8324221220657863e-7;
	a[7] = -6.7778497812918602e-6;

	double re_zn = 1.0;
	double im_zn = 0.0;

	for (i=0; i<8; i++)
	{
		// a[i] += 0.5 / ((double) (i+1));
		rep += a[i] * re_zn;
		imp += a[i] * im_zn;
		
		tmp = re_q * re_zn - im_q * im_zn;
		im_zn = re_q * im_zn + im_q * re_zn;
		re_zn = tmp;
	}
	
	*prep = rep;
	*pimp = imp;
}
			  
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

	/* z * ln gamma (1/(1-z)) */
	rep += re_q * re_lng - im_q * im_lng;
	imp += re_q * im_lng + im_q * re_lng;

	/* log (1-z) */
	gsl_sf_result m_lnz, t_lnz;
	gsl_sf_complex_log_e (re_one_minus_z, im_one_minus_z, &m_lnz, &t_lnz);
	double re_ln_omz = m_lnz.val;
	double im_ln_omz = t_lnz.val;
	
	/* (1/z) * (1/(1-z)) */
	double re_zz = re_one_over_z * re_one_over_one_minus_z - im_one_over_z * im_one_over_one_minus_z;
	double im_zz = re_one_over_z * im_one_over_one_minus_z + im_one_over_z * re_one_over_one_minus_z;

	/* (1+z) / (z(1-z)) */
	double re_tz = (1.0+re_q) * re_zz - im_q * im_zz;
	double im_tz = (1.0+re_q) * im_zz + im_q * re_zz;

	/* add 0.5 * log (1-z) * (1+z)/(z(1-z)) */
	rep += 0.5 * (re_ln_omz * re_tz - im_ln_omz * im_tz);
	imp += 0.5 * (re_ln_omz * im_tz + im_ln_omz * re_tz);
	
	/* add 1/(1-z) */
	rep += re_one_over_one_minus_z;
	imp += im_one_over_one_minus_z;

	rep += 0.577215664901532 * re_one_over_z;
	imp += 0.577215664901532 * im_one_over_z;
	
	*prep = rep;
	*pimp = imp;
}

static double alpha_series (double re_q, double im_q)
{
	double rep, imp;
	
#if AT_INFINITY
	tmp = 1.0 / (re_q *re_q + im_q * im_q);
	reo = re_q * tmp;
	imo = -im_q * tmp;
	re_q = reo;
	im_q = imo;
#endif
	
	alpha_series_c (re_q, im_q, &rep, &imp);

	// const_series_c (re_q, im_q, &rep, &imp);
	
#if 0
	double reo = rep * re_q - imp * im_q;
	double imo = rep * im_q + imp * re_q;

	rep = reo;
	imp = imo;
#endif

	return sqrt (rep*rep+imp*imp);
	// return log (sqrt (rep*rep+imp*imp));
	// return imp;
	// return (atan2 (imp, rep) + M_PI) / (2.0*M_PI);
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
