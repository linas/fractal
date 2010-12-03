/*
 * swap.C
 *
 * FUNCTION:
 * Integral the permuation group of continued fractions
 * Expect to get Riemann zeta in the Gauss map case
 * and that is what we seem to get ... need high integration 
 * order though to get anything on the r=1/2 axis ... 
 *
 * Results written up in yarh.lyx
 *
 * Linas Feb 2005
 */ 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

inline long double swap12 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	if (1.0e-10 > r1) return r1;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	r2 += a1;
	r1 = 1.0L / r2;
	r1 += a2;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap13 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	if (1.0e-10 > r1) return 0.69314718;  // "average" value
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-10 > r2) return r2;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;

	r3 += a1;

	r2 = 1.0L / r3;
	r2 += a2;
	
	r1 = 1.0L / r2;
	r1 += a3;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap23 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	if (1.0e-10 > r1) return 0.0;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-10 > r2) return 0.0;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;

	r3 += a2;

	r2 = 1.0L / r3;
	r2 += a3;
	
	r1 = 1.0L / r2;
	r1 += a1;
	
	x = 1.0L / r1;
	return x;
}

inline long double swap14 (long double x)
{
	long double ox = 1.0L/x;
	long double a1 = floorl(ox);
	long double r1 = ox - a1;
	if (1.0e-10 > r1) return 0.0;
	ox = 1.0L/r1;
	long double a2 = floorl(ox);
	long double r2 = ox - a2;
	if (1.0e-10 > r2) return 0.0;
	ox = 1.0L/r2;
	long double a3 = floorl(ox);
	long double r3 = ox - a3;
	if (1.0e-10 > r3) return 0.0;
	ox = 1.0L/r3;
	long double a4 = floorl(ox);
	long double r4 = ox - a4;

	r4 += a1;

	r3 = 1.0L / r4;
	r3 += a3;
	
	r2 = 1.0L / r3;
	r2 += a2;
	
	r1 = 1.0L / r2;
	r1 += a4;
	
	x = 1.0L / r1;
	return x;
}

/* The integrand, which is swap(x) * x^s */
void grand (long double x, long double sre, long double sim, 
					 long double *pre, long double *pim)
{
#if 0
	// This is the basic case, which gives Riemann exactly
	long double ox = 1.0L/x;
	long double sw = ox - floorl(ox);
#else
	long double sw = swap12 (x);
	// long double sw = swap13 (x);
	// long double sw = swap23 (x);

	// This is used for the sanity-check, "shadow" graph
	// long double sw = x;
#endif

	long double lnx = logl (x);
	long double ire = expl (sre*lnx);
	// long double ire = 1.0L/ sqrtl (x);
	long double iim = ire;
	long double phi = sim*lnx;
	ire *= cosl (phi);
	iim *= sinl (phi);
	
	ire *= sw;
	iim *= sw;

	*pre = ire;
	*pim = iim;
}


/* compute single integral of the integrand.
 * actually compute 
 * zeta = s/(s-1) - s \int_0^1 swap(x) x^{s-1} dx 
 */
void gral(int nsteps, long double sre, long double sim, 
					 long double *pre, long double *pim)
{
	int i;

	long double step = 1.0L / ((long double) nsteps);
	
	long double sum_re= 0.0L;
	long double sum_im= 0.0L;
	long double x = 1.0L - 0.5*step;

	long double r = RAND_MAX;
	r = 1.0L / r;

	/* integrate in a simple fashion */
	int nh = 0;
	for (i=0; i<nsteps; i++)
	{
// #define DO_RAND
#ifdef DO_RAND
		int nr = rand();
		if (0 == nr) continue;
		x = (long double) nr;
		x *= r;
#endif
		long double val_re, val_im;
		grand (x, sre-1.0, sim, &val_re, &val_im);
		x -= step;
		sum_re += val_re;
		sum_im += val_im;
		nh ++;
	}

	/* Divide by the actual number of samples */
	step = 1.0L / ((long double) nh);
	sum_re *= step;
	sum_im *= step;

#if 0
	// the trivial case
	sum_re = sre+1.0;
	sum_im = sim;
	long double t = sum_re*sum_re + sum_im*sum_im;
	sum_re /= t;
	sum_im = -sum_im / t;
#endif

	/* mult by s */
	long double tmp = sum_re * sre - sum_im*sim;
	sum_im = sum_im *sre + sum_re *sim;
	sum_re = tmp;

	/* compute 1/(s-1) */
	sre -= 1.0L;
	long double v = sre*sre+sim*sim;
	sre /= v;
	sim = -sim/v;

	/* subtract */
	sum_re = sre - sum_re;
	sum_im = sim - sum_im;

	/* s/(s-1) = 1/(s-1) + 1 so add 1 now */
	sum_re += 1.0L;

	*pre = sum_re;
	*pim = sum_im;
}

double
rswap (long double sre, long double sim, int itermax)
{
	long double zre, zim;

	gral (itermax, sre, sim, &zre, &zim);
	long double mag = zre*zre+zim*zim;
	mag = sqrt (mag);
	// return mag;
	long double  phase = atan2l(zim, zre);
	phase /= 2.0L * M_PI;
	phase += 0.5L;
	return phase;
}

// DECL_MAKE_HEIGHT(rswap)
/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
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
   double im_end = im_center - width * ((double) sizey) / (2.0 * (double) sizex);

	printf ("re=(%g,%g)\n", re_start, re_start+width);
	printf ("im=(%g,%g)\n", im_end, im_start);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) { printf(" start row %d\n", i); fflush (stdout); }
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = rswap (re_position, im_position, itermax);
         glob [i*sizex +j] = phi;

			// draw vertical lines showing crit strip 
			if ((re_position <= 0.0) && (0.0<re_position+delta))
         	glob [i*sizex +j] = -1.0;
			if ((re_position <= 0.5) && (0.5<re_position+delta))
         	glob [i*sizex +j] = -1.0;
			if ((re_position <= 1.0) && (1.0<re_position+delta))
         	glob [i*sizex +j] = -1.0;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
