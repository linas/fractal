/*
 * cutoff.C
 *
 * FUNCTION:
 * Explore spectral analysis of the interior of the Mandelbrot set
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

#include "brat.h"

/*-------------------------------------------------------------------*/
/* This routine does a spectral analysis for the Mandelbrot set iterator.
 * That is, it computes a reimann-zeta-like sum of things like the modulus
 * (dirichlet series, to be precise)
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
   double	re, im, tmp, mod, modulus=0.0;
   double	dre, dim, dmod;
   double	ddre, ddim, ddmod;
   double	zpre, zpim, zppre, zppim;
   double	mp, mpp;
   int		loop;
   double 	omod=0.0,  theta;
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
   rp        = (double *) malloc ((itermax+1)*sizeof (double));
   rpp       = (double *) malloc ((itermax+1)*sizeof (double));
   rppp      = (double *) malloc ((itermax+1)*sizeof (double));
   rpppp     = (double *) malloc ((itermax+1)*sizeof (double));

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
   printf ("itermax=%d tau=%g 1/tau=%g sum_n=%g sum_np=%g sum_npp=%g\n", 
            itermax, tau, 1.0/tau, sum_n, sum_np, sum_npp);
   printf (" n^2=%g 2n^3=%g\n", sum_n*sum_n, 2.0*sum_n*sum_n*sum_n);

   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0 / log(2.0);
   

   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   

// sizey=1;
// im_start=0.0;
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
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
         dre = 1.0;
         dim = 0.0;
         dmod = 0.0;
         ddre = 0.0;
         ddim = 0.0;
         ddmod = 0.0;
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
            zppre = re*ddre + im*ddim;   /* divergence */
            zppim = re*ddim - im*ddre;   /* curl */
            zppre *= omod;
            zppim *= omod;

            /* compute zprime/z */
            zpre = re*dre + im*dim;   /* divergence */
            zpim = re*dim - im*dre;   /* curl */
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
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

#if 0
         modulus = sqrt (modulus);
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         frac = log (log (modulus)) * tl;

         glob [i*sizex +j] = modulus / sum_n;
         glob [i*sizex +j] = sum_mod / sum_n;
         glob [i*sizex +j] = sum_mod - modulus;

          /* the following computes an almost-flat, divergence free thing */
         glob [i*sizex +j] = modulus / sum_n;
         glob [i*sizex +j] /= (sqrt (re_position*re_position+im_position*im_position));

         /* the interesting one is the z-prime-prime */
         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         glob [i*sizex +j] = modulus / sum_n;
         // glob [i*sizex +j] -= 0.25 * exp( -0.75 * log((re_position-0.25)*(re_position-0.25)+im_position*im_position));

         /* --------------------------------------------------------- */
         /* the interesting one is the z-prime-prime */
         /* here we use a taylor expansion to extrapolate to tau=0 */
         /* first, we need the derivatives of modulus w.r.t tau */
         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         mp = (sum_ddrep * sum_ddre + sum_ddimp * sum_ddim) / modulus;
         mpp  = sum_ddrep * sum_ddrep + sum_ddre * sum_ddrepp;
         mpp += sum_ddimp * sum_ddimp + sum_ddim * sum_ddimpp - mp*mp;
         mpp /= modulus;
         
         /* next, we need derivatives of m/n w.r.t tau */
         mod = modulus / sum_n;
         dmod = (mp - mod * sum_np) / sum_n;
         ddmod = (mpp - 2.0 * dmod * sum_np - mod * sum_npp) / sum_n;

         /* finally the taylor expansion */
         glob [i*sizex +j] = mod - tau* (dmod - 0.5 * tau * ddmod);

         // glob [i*sizex +j] -= 0.25 * exp( -0.75 * log((re_position-0.25)*(re_position-0.25)+im_position*im_position));

         /* --------------------------------------------------------- */
         /* OK, lets do the taylor expansion for just-plain Z */
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         mp = (sum_rep * sum_re + sum_imp * sum_im) / modulus;
         mpp  = sum_rep * sum_rep + sum_re * sum_repp;
         mpp += sum_imp * sum_imp + sum_im * sum_impp - mp*mp;
         mpp /= modulus;
         
         /* next, we need derivatives of m/n w.r.t tau */
         mod = modulus / sum_n;
         dmod = (mp - mod * sum_np) / sum_n;
         ddmod = (mpp - 2.0 * dmod * sum_np - mod * sum_npp) / sum_n;

         /* finally the taylor expansion */
         glob [i*sizex +j] = mod - tau* (dmod - 0.5 * tau * ddmod);
// printf ("%9.6g	%9.6g	%9.6g\n", re_position, glob[i*sizex+j], mod);

         /* ok, this part should be the divergent part ... */
         theta = 0.5 * atan2 (-im_position, 0.25-re_position);
         mod = (re_position-0.25)*(re_position-0.25)+im_position*im_position;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += 0.5;
         tmp = sqrt(re*re+im*im);
         if (0.5 < tmp) tmp = 0.5;

         glob [i*sizex +j] -= tmp;

         /* --------------------------------------------------------- */
         /* the interesting one is the z-prime-prime */
         /* here we use a taylor expansion to extrapolate to tau=0 */
         /* first, we need the drerivatives of modulus w.r.t tau */
         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         mp = (sum_ddrep * sum_ddre + sum_ddimp * sum_ddim) / modulus;
         mpp  = sum_ddrep * sum_ddrep + sum_ddre * sum_ddrepp;
         mpp += sum_ddimp * sum_ddimp + sum_ddim * sum_ddimpp - mp*mp;
         mpp /= modulus;
         
         /* finally the taylor expansion */
         /* subtract the main-body divergence */
         tmp = 0.25 * exp( -0.75 * log((re_position-0.25)*(re_position-0.25)+im_position*im_position));
         
         glob [i*sizex +j] = (modulus-tmp*sum_n);
         glob [i*sizex +j] -= tau* ((mp-tmp*sum_np) - 0.5 * tau * (mpp-tmp*sum_npp));
#endif

         /* --------------------------------------------------------- */
         /* OK, lets do the taylor expansion for just-plain Z */
         /* here, we subtract the leading divergence 
          */
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         mp = (sum_rep * sum_re + sum_imp * sum_im) / modulus;
         mpp  = sum_rep * sum_rep + sum_re * sum_repp;
         mpp += sum_imp * sum_imp + sum_im * sum_impp - mp*mp;
         mpp /= modulus;
         
         /* finally the taylor expansion */
         glob [i*sizex +j] = modulus - tau* (mp - 0.5 * tau * mpp);

         theta = 0.5 * atan2 (-im_position, 0.25-re_position);
         mod = (re_position-0.25)*(re_position-0.25)+im_position*im_position;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += 0.5;
         tmp = sqrt(re*re+im*im);
         if (0.5 < tmp) tmp = 0.5;

         glob [i*sizex +j] -= tmp * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));

#if 0
         /* --------------------------------------------------------- */
         /* the interesting one is the z-prime-prime */
         /* here we use a taylor expansion to extrapolate to tau=0 */
         /* first, we need the derivatives of modulus w.r.t tau */
         /* we compute the phase */
         phi = atan2 (sum_ddim, sum_ddre);
         modulus = 1.0 / sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         phip = (sum_ddre * sum_ddimp - sum_ddrep * sum_ddim) * modulus;
         phipp = (sum_ddre * sum_ddimp - sum_ddrep * sum_ddim) * modulus;
         phipp *=  -2.0* (sum_ddre * sum_ddrep + sum_ddim * sum_ddimp) * modulus;
         phipp += (sum_ddre * sum_ddimpp - sum_ddrepp * sum_ddim) * modulus;

         glob [i*sizex +j] = (phi + M_PI)/(2.0*M_PI);
         // glob [i*sizex +j] = (phi - tau* (phip - 0.5 * tau * phipp) +M_PI)/(2.0*M_PI);

         theta = -1.5 * atan2 (-im_position, 0.25-re_position);
         mod = (re_position-0.25)*(re_position-0.25)+im_position*im_position;
         mod = pow (mod, -0.75);
         re = 0.25 * mod * cos(theta);
         im = 0.25 * mod * sin(theta);

         glob [i*sizex +j] =  (sum_ddre/sum_n-re)*(sum_ddre/sum_n-re);
         glob [i*sizex +j] += (sum_ddim/sum_n-im)*(sum_ddim/sum_n-im);
         glob [i*sizex +j] = sqrt (glob[i*sizex+j]);

         /* --------------------------------------------------------- */
         theta = 0.5 * atan2 (-im_position, 0.25-re_position);
         mod = (re_position-0.25)*(re_position-0.25)+im_position*im_position;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += re_position;
         im += im_position;

         glob [i*sizex +j] =  sqrt (re*re +im*im);
#endif
         /* --------------------------------------------------------- */
         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
