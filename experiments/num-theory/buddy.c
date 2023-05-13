/*
 * buddy.c
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

/*-------------------------------------------------------------------*/
/* This routine does a spectral analysis for the Mandelbrot set iterator.
 * That is, it computes a reimann-zeta-like sum of things like the modulus
 * (dirichlet series, to be precise)
 */

void
bud_sum (int sizex, int itermax)
{
	int		i,j;
	long double	re_start, delta;
	long double	re_position, im_position;
	long double	re, im, tmp, modulus=0.0;
	long double	dre, dim, dmod;
	long double	ddre, ddim, ddmod;
	long double	zpre, zpim, zppre, zppim;
	int		loop;
	long double 	omod=0.0;
	long double 	escape_radius = 1.0e30;
	long double 	ren, tl;
	long double	tau;
	long double	*regulator, *rp, *rpp, *rppp, *rpppp;
	long double	sum_n, sum_np, sum_npp, sum_nppp, sum_npppp;
	long double	sum_re, sum_im, sum_mod;
	long double	sum_rep, sum_imp, sum_modp;
	long double	sum_repp, sum_impp, sum_modpp;
	long double	sum_dre, sum_dim, sum_dmod;
	long double	sum_ddre, sum_ddim, sum_ddmod;
	long double	sum_ddrep, sum_ddimp, sum_ddmodp;
	long double	sum_ddrepp, sum_ddimpp, sum_ddmodpp;
	long double	sum_zpre, sum_zpim, sum_zpmod;
	long double	sum_zppre, sum_zppim, sum_zppmod;

	/* first, compute the regulator, so that the itermax'th iteration makes 
	 * a negligable contribution (about 1e-30) */
	tau = 16.0 / ((long double) itermax);

	/* set up smooth ramp 
	 * regulator is exponential, and rp is derivative w.r.t. tau,
	 * rpp is second deriv. w.r.t. tau */
	sum_n = sum_np = sum_npp = sum_nppp = sum_npppp = 0.0;
	regulator = (long double *) malloc ((itermax+1)*sizeof (long double));
	rp		  = (long double *) malloc ((itermax+1)*sizeof (long double));
	rpp		 = (long double *) malloc ((itermax+1)*sizeof (long double));
	rppp		= (long double *) malloc ((itermax+1)*sizeof (long double));
	rpppp	  = (long double *) malloc ((itermax+1)*sizeof (long double));

	for (i=0; i<itermax; i++) 
	{
		tmp = - (long double) (i*i);
		regulator[i] = exp (tmp * tau*tau);
		rp[i] = 2.0 * tau * tmp * regulator[i];
		rpp[i] = 2.0 * tmp * (regulator[i] + tau * rp[i]);
		sum_n += regulator[i];
		sum_np += rp[i];
		sum_npp += rpp[i];
	}
	printf ("#\n# interior of mandelbrot bud\n#\n");
	printf ("# itermax=%d tau=%llg 1/tau=%llg sum_n=%llg tau*sum_n=%llg\n", 
				itermax, tau, 1.0/tau, sum_n, tau*sum_n);
	printf ("# sum_np=%g sum_npp=%llg\n", sum_np, sum_npp);
	printf ("#  n^2=%llg 2n^3=%llg\n", sum_n*sum_n, 2.0*sum_n*sum_n*sum_n);
	printf ("#  n - tau* (np - 0.5 * tau * npp) = %llg\n",  
	         sum_n - tau* (sum_np - 0.5 * tau * sum_npp));
	printf ("#  tau*(n - tau* (np - 0.5 * tau * npp)) = %llg\n",  
	         tau*(sum_n - tau* (sum_np - 0.5 * tau * sum_npp)));

	ren = log( log (escape_radius)) / log(2.0);
	tl = 1.0 / log(2.0);

	/* start at the center of the bud which is at -1.0 */
	delta = 0.25 / (long double) sizex;
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
		long double re_c = re_position;
		long double im_c = im_position;
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

		long double q = 4.0L*(re_position+1.0L);
		modulus -= 3.0 - 8.96e-06;
		modulus -= 7.5 * q;
		modulus -= 10.5*q*q;
		modulus -= 20.5*q*q*q;
		modulus -= 0.0*q*q*q*q;
		modulus -= 65.0*q*q*q*q*q;

		if (j != 0) modulus /= q*q*q*q*q*q;

		sum_ddre -= 3.0;
		sum_ddre -= 7.5*q;
		sum_ddre -= 10.8*q*q;
		sum_ddre -= 18.0*q*q*q;
		sum_ddre -= 13.0*q*q*q*q;

		if (j != 0) sum_ddre /= q*q*q*q*q;

		printf ("%d	%llg	%llg	%llg\n", 
			j, re_position, q, sum_ddim);

		/* --------------------------------------------------------- */
		re_position += delta;
	}
}

int main (int argc, char * argv[])
{
	if (3 > argc)
	{
		fprintf (stderr, "Usage: %s <npoints> <itermax>\n", argv[0]);
		exit (1);
	}

	int npts = atoi (argv[1]);
	int itermax = atoi (argv[2]);

	bud_sum (npts, itermax);
}

/* --------------------------- END OF LIFE ------------------------- */
