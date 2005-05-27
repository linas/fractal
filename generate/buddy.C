/*
 * buddy.C
 *
 * FUNCTION:
 * Explore spectral analysis of the interior of the 
 * western bud of the Mandelbrot set
 *
 * HISTORY:
 * Based on version from 2000
 * more stuff -- October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coord-xforms.h"

/*-------------------------------------------------------------------*/
/* This routine does a spectral analysis for the Mandelbrot set iterator.
 * That is, it computes a reimann-zeta-like sum of things like the modulus
 * (dirichlet series, to be precise)
 */

void
bud_sum (int sizex, int itermax)
{
	int		i,j;
	double	re_start, delta;
	double	re_position, im_position;
	double	re, im, tmp, modulus=0.0;
	double	dre, dim, dmod;
	double	ddre, ddim, ddmod;
	double	zpre, zpim, zppre, zppim;
	int		loop;
	double 	omod=0.0;
	double 	escape_radius = 1.0e30;
	double 	ren, tl;
	double	tau;
	double	*regulator, *rp, *rpp, *rppp, *rpppp;
	double	sum_n, sum_np, sum_npp, sum_nppp, sum_npppp;
	double	sum_re, sum_im, sum_mod;
	double	sum_rep, sum_imp, sum_modp;
	double	sum_repp, sum_impp, sum_modpp;
	double	sum_dre, sum_dim, sum_dmod;
	double	sum_ddre, sum_ddim, sum_ddmod;
	double	sum_ddrep, sum_ddimp, sum_ddmodp;
	double	sum_ddrepp, sum_ddimpp, sum_ddmodpp;
	double	sum_zpre, sum_zpim, sum_zpmod;
	double	sum_zppre, sum_zppim, sum_zppmod;

	/* first, compute the regulator, so that the itermax'th iteration makes 
	 * a negligable contribution (about 1e-30) */
	tau = 16.0 / ((double) itermax);

	/* set up smooth ramp 
	 * regulator is exponential, and rp is derivative w.r.t. tau,
	 * rpp is second deriv. w.r.t. tau */
	sum_n = sum_np = sum_npp = sum_nppp = sum_npppp = 0.0;
	regulator = (double *) malloc ((itermax+1)*sizeof (double));
	rp		  = (double *) malloc ((itermax+1)*sizeof (double));
	rpp		 = (double *) malloc ((itermax+1)*sizeof (double));
	rppp		= (double *) malloc ((itermax+1)*sizeof (double));
	rpppp	  = (double *) malloc ((itermax+1)*sizeof (double));

	for (i=0; i<itermax; i++) 
	{
		tmp = - (double) (i*i);
		regulator[i] = exp (tmp * tau*tau);
		rp[i] = 2.0 * tau * tmp * regulator[i];
		rpp[i] = 2.0 * tmp * (regulator[i] + tau * rp[i]);
		sum_n += regulator[i];
		sum_np += rp[i];
		sum_npp += rpp[i];
	}
	printf ("# itermax=%d tau=%g 1/tau=%g sum_n=%g tau*sum_n=%g\n", 
				itermax, tau, 1.0/tau, sum_n, tau*sum_n);
	printf ("# sum_np=%g sum_npp=%g\n", sum_np, sum_npp);
	printf ("#  n^2=%g 2n^3=%g\n", sum_n*sum_n, 2.0*sum_n*sum_n*sum_n);
	printf ("#  n - tau* (np - 0.5 * tau * npp) = %g\n",  
	         sum_n - tau* (sum_np - 0.5 * tau * sum_npp));
	printf ("#  tau*(n - tau* (np - 0.5 * tau * npp)) = %g\n",  
	         tau*(sum_n - tau* (sum_np - 0.5 * tau * sum_npp)));

	ren = log( log (escape_radius)) / log(2.0);
	tl = 1.0 / log(2.0);

	/* start at the center of the bud which is at -1.0 */
	delta = 0.25 / (double) sizex;
	re_start = -1.0;
	
	im_position = 0.0;
	re_position = re_start;
	for (j=0; j<sizex; j++)
	{
		sum_re = sum_im = sum_mod = 0.0;
		sum_rep = sum_imp = sum_modp = 0.0;
		sum_repp = sum_impp = sum_modpp = 0.0;
		sum_dre = sum_dim = sum_dmod = 0.0;
		sum_ddre = sum_ddim = sum_ddmod = 0.0;
		sum_ddrep = sum_ddimp = sum_ddmodp = 0.0;
		sum_ddrepp = sum_ddimpp = sum_ddmodpp = 0.0;
		sum_zpre = sum_zpim = sum_zpmod = 0.0;
		sum_zppre = sum_zppim = sum_zppmod = 0.0;
		re = re_position;
		im = im_position;
		double re_c = re_position;
		double im_c = im_position;
		dre = 1.0;
		dim = 0.0;
		dmod = 0.0;
		ddre = 0.0;
		ddim = 0.0;
		ddmod = 0.0;

		re = re_c;
		im = im_c;
		modulus = (re*re + im*im);
		for (loop=1; loop <itermax; loop++) 
		{
			sum_re += re * regulator [loop];
			sum_im += im * regulator [loop];
			sum_rep += re * rp [loop];
			sum_imp += im * rp [loop];
			sum_repp += re * rpp [loop];
			sum_impp += im * rpp [loop];
			// sum_mod += sqrt(modulus) * regulator [loop];

			/* sum over first derivative z-prime */
			sum_dre += dre * regulator [loop];
			sum_dim += dim * regulator [loop];
			// sum_dmod += sqrt(dmod) * regulator [loop];

			/* sum over second derivative z-prime-prime*/
			sum_ddre += ddre * regulator [loop];
			sum_ddim += ddim * regulator [loop];
			sum_ddrep += ddre * rp [loop];
			sum_ddimp += ddim * rp [loop];
			sum_ddrepp += ddre * rpp [loop];
			sum_ddimpp += ddim * rpp [loop];
			// sum_ddmod += sqrt(ddmod) * regulator [loop];

			omod = 1.0 / modulus;

			/* compute zprimeprime/z */
			zppre = re*ddre + im*ddim;	/* divergence */
			zppim = re*ddim - im*ddre;	/* curl */
			zppre *= omod;
			zppim *= omod;

			/* compute zprime/z */
			zpre = re*dre + im*dim;	/* divergence */
			zpim = re*dim - im*dre;	/* curl */
			zpre *= omod;
			zpim *= omod;

			/* sum over first and second derivative z-prime-prime / z*/
			sum_zpre += zpre * regulator [loop];
			sum_zpim += zpim * regulator [loop];
			sum_zppre += zppre * regulator [loop];
			sum_zppim += zppim * regulator [loop];

			/* compute second derivative */
			tmp = 2.0 * (re*ddre - im*ddim + dre*dre - dim*dim);
			ddim = 2.0 * (re*ddim + im*ddre + 2.0 * dre*dim);
			ddre = tmp;
			ddmod = ddre*ddre+ddim*ddim;

			/* compute infinitessimal flow */
			tmp = 2.0 * (re*dre - im*dim) +1.0;
			dim = 2.0 * (re*dim + im*dre);
			dre = tmp;
			dmod = dre*dre+dim*dim;

			/* compute iterate */
			tmp = re*re - im*im + re_c;
			im = 2.0*re*im + im_c;
			re = tmp;
			modulus = (re*re + im*im);
			if (modulus > escape_radius*escape_radius) break;
		}	 

		/* --------------------------------------------------------- */
		/* The interesting one is the z-prime-prime. In the main 
		 * bud to the west, its finite. */
		modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);

		printf ("%d	%g	%g\n", j, re_position, modulus);

		/* --------------------------------------------------------- */
		re_position += delta;
	}
}

int main (int argc, char * argv[])
{
	bud_sum (20, 20000);
}

/* --------------------------- END OF LIFE ------------------------- */
