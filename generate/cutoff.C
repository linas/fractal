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
#include "coord-xforms.h"

/* return real part of mobius x-form on the poincare disk */
static inline double 
remob (double a, double b, double c, double d, double x, double y)
{
	double nr = (a+d)*x + (c-b)*y +b+c;
	double ni = (b-c)*x + (a+d)*y +a-d;
	double dr = (b+c)*x + (a-d)*y +a+d;
	double di = (d-a)*x + (b+c)*y +c-b;

	double rem = (nr*dr + ni*di) / (dr*dr+di*di);
	return rem;
}
              
/* return imaginary part of mobius x-form on the poincare disk */
static inline double 
immob (double a, double b, double c, double d, double x, double y)
{
	double nr = (a+d)*x + (c-b)*y +b+c;
	double ni = (b-c)*x + (a+d)*y +a-d;
	double dr = (b+c)*x + (a-d)*y +a+d;
	double di = (d-a)*x + (b+c)*y +c-b;

	double imm = (ni*dr - nr*di) / (dr*dr+di*di);
	return imm;
}
              
/* The following is the set of defines needed to be enabled to
 * get the flattening of the main cardioid into q-disk coords to work.
 */
// #define MAIN_CARDIOID_MODULAR_FORM
#ifdef MAIN_CARDIOID_MODULAR_FORM
// #define Q_SERIES_MOBIUS
#define CIRCLE_COORDS
#define FLATTEN_CARDIOID_MAP
// #define ZPP_MODULUS_DIVERGENCE_FREE
#define COMPLEX_ZPP_MINUS_DIVERGENCE
#endif

/* The following sets up the complex divisor-sum-like thingy
 * in the main bud. 
 */
#define BUD_COMPLEX_FORM
#ifdef BUD_COMPLEX_FORM
#define Q_SERIES_MOBIUS
#define CIRCLE_TO_BUD
#define PLAIN_Z_PRIME_PRIME
#endif

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
	printf ("#\n#interiopr of mandelbrot\n#\n");
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
			double re_c = re_position;
			double im_c = im_position;
         dre = 1.0;
         dim = 0.0;
         dmod = 0.0;
         ddre = 0.0;
         ddim = 0.0;
         ddmod = 0.0;

// #define Q_SERIES_MOBIUS
#ifdef Q_SERIES_MOBIUS
			/* First, make a map from q-series coords to the 
			 * upper half-plane, then apply the mobius x-form, 
			 * and then go back to the q-series coords */
			double tau_re, tau_im;

#ifdef ROT
			poincare_disk_to_plane_coords (re_c, im_c, &tau_re, &tau_im);
			double a=1, b=1, c=0, d=1;
			mobius_xform (a,b,c,d, tau_re, tau_im, &tau_re, &tau_im);
			plane_to_poincare_disk_coords (tau_re, tau_im, &re_c, &im_c);
#endif

			// poincare_disk_to_plane_coords (-im_c, re_c, &tau_re, &tau_im);
			// poincare_disk_to_plane_coords (re_c, im_c, &tau_re, &tau_im);
			q_disk_to_plane_coords (re_c, im_c, &tau_re, &tau_im);
			mobius_xform (2,0,0,1, tau_re, tau_im, &tau_re, &tau_im);

			double a=1, b=0, c=1, d=1;
			// double a=3, b=1, c=0, d=1;
			// double a=0, b=-1, c=1, d=0;
			mobius_xform (a,b,c,d, tau_re, tau_im, &tau_re, &tau_im);

			/* compute the modular scaling factor */
c=0;
			double re_var = c*tau_re+d;
			double im_var = c*tau_im;
			double var = re_var*re_var + im_var*im_var;

			double re_2v = re_var*re_var - im_var*im_var;
			double im_2v = 2.0 * re_var* im_var;

			double re_iv = re_var / var;
			double im_iv = -im_var / var;

			double angle = atan2 (im_var, re_var);
			// angle += 0.5*M_PI;
			double sca = 1.0;
			double re_sca = pow (var, 0.5*sca) * cos (sca*angle);
			double im_sca = pow (var, 0.5*sca) * sin (sca*angle);

			plane_to_q_disk_coords (tau_re, tau_im, &re_c, &im_c);
			// plane_to_poincare_disk_coords (tau_re, tau_im, &re_c, &im_c);

#endif /* Q_SERIES_MOBIUS */

// #define CIRCLE_TO_BUD
#ifdef CIRCLE_TO_BUD
			/* Map a unit disk centered at origin to the main bud
			 * located to the west */
			re_c *= 0.25;
			im_c *= 0.25;
			re_c -= 1.0;
#endif /* CIRCLE_TO_BUD */

#ifdef POINCARE_CIRCLE_MOBIUS
			/* This is the mobius map for the poincare disk, which is 
			 * incorrect for the punctured disk aka q-series disk */
			double a,b,c,d;
			a = 1; b=1; c=0; d=1;
			double xp = remob (a,b,c,d, re_c, im_c);
			double yp = immob (a,b,c,d, re_c, im_c);
			re_c = xp;
			im_c = yp;
#endif /* MOBIUS */

//xx #define CIRCLE_COORDS
#ifdef CIRCLE_COORDS
			/* This map takes a circle and maps it to a rectangle,
			 * i.e. radius-theta coords, needed as input to cardioid */

			/* A left-to-right inversion is needed to get the image in the
			 * same coordinate frame as the q-series coords */
			re_c = -re_c;

			double rr = sqrt(re_c*re_c + im_c*im_c);
			double theta = atan2 (im_c, re_c);
			theta /= M_PI;
			re_c = theta;
			im_c = rr;
#endif /* CIRCLE_COORDS */

//xx #define FLATTEN_CARDIOID_MAP
#ifdef FLATTEN_CARDIOID_MAP
         /* Map to cardiod lam(1-lam) 
			 * Input to this thing is assumed to be a ractangle, 
			 * going from x= -1.0 to 1.0 and y= 0 to 1
			 * which gets mapped to cardiod with y=1 at the edge,
			 * and x=0 at the left side of cardioid
			 */
         double r = -im_c;
         double phi = M_PI*re_c;

         /* works pretty well
         r *= r;
         r = 1.0 + (r-1.0)*sin(0.5*phi)*sin(0.5*phi);
         */
         // r -= 1.0;
         // r *= sin(0.5*phi)*sin(0.5*phi);
         // r *= (0.5*phi)*sin(0.5*phi);
         // r *= (0.5*phi)* (0.5*phi);
        
         // r += 1.0;
         re_c = 0.5 * r * (cos (phi) - 0.5 * r * cos (2.0*phi));
         im_c = 0.5 * r * (sin (phi) - 0.5 * r * sin (2.0*phi));
#endif /* FLATTEN_CARDIOID_MAP */

#ifdef ALT_FLATTEN_DOESNT_WORK_AT_ALL
         double r = im_position;
         double phi = M_PI*re_position;

			r *= 1.0- cos(phi);

			re_c = 0.25+ r*cos(phi);
			im_c = r*sin(phi);

#endif /* ALT_FLATTEN */

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
            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

#if MISC_SIMPLE_THINGS
         modulus = sqrt (modulus);
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         frac = log (log (modulus)) * tl;

         glob [i*sizex +j] = modulus / sum_n;
         glob [i*sizex +j] = sum_mod / sum_n;
         glob [i*sizex +j] = sum_mod - modulus;

          /* the following computes an almost-flat, divergence free thing */
         glob [i*sizex +j] = modulus / sum_n;
         glob [i*sizex +j] /= (sqrt (re_c*re_c+im_c*im_c));
#endif

// #define PLAIN_Z_PRIME_PRIME
#ifdef PLAIN_Z_PRIME_PRIME
         /* --------------------------------------------------------- */
         /* The interesting one is the z-prime-prime. In the main 
			 * bud to the west, its finite. */

#ifdef LEADING_TERMS
         // modulus -= 3.0;
         // modulus -= 7.5*q;
         // modulus -= 10.5*q*q;
         // modulus -= 20.5*q*q*q;
         // modulus -= 0.0*q*q*q*q;
         // modulus -= 65.0*q*q*q*q*q;

         sum_ddre -= 3.0;
			double qre = 4.0*(re_c+1.0);
			double qim = 4.0*im_c;
         sum_ddre -= 7.5 * qre;
         sum_ddim -= 7.5 * qim;

			double q2re = qre*qre - qim*qim;
			double q2im = 2.0*qre*qim;

         sum_ddre -= 10.8 * q2re;
         sum_ddim -= 10.8 * q2im;

			double q3re = q2re*qre - q2im*qim;
			double q3im = q2re*qim + qre*q2im;

         sum_ddre -= 19.0 * q3re;
         sum_ddim -= 19.0 * q3im;
#endif /* LEADING_TERMS */

         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);

// #define MODULAR_FORM_CORRECTIONS
#ifdef MODULAR_FORM_CORRECTIONS
         glob [i*sizex +j] = sum_ddim;
         glob [i*sizex +j] = modulus;
         glob [i*sizex +j] = sum_ddre;

         glob [i*sizex +j] = sum_ddre*re_2v - sum_ddim*im_2v;
         glob [i*sizex +j] = sum_ddre*re_var - sum_ddim*im_var;
         glob [i*sizex +j] = sum_ddre*re_iv - sum_ddim*im_iv;
         glob [i*sizex +j] = sum_ddre*re_sca - sum_ddim*im_sca;

         double rem = sum_ddre*re_sca - sum_ddim*im_sca;
         double imm = sum_ddre*im_sca + sum_ddim*re_sca;
         modulus = sqrt (rem*rem + imm*imm);
#endif /* MODULAR_FORM_CORRECTIONS */

         glob [i*sizex +j] = fabs(sum_ddim);
         glob [i*sizex +j] = sum_ddre;
         glob [i*sizex +j] = modulus;
#endif

#ifdef PLAIN_Z_PRIME_PRIME_NORMALIZED
         /* --------------------------------------------------------- */
         /* the interesting one is the z-prime-prime */
         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         glob [i*sizex +j] = modulus / sum_n;
         // glob [i*sizex +j] -= 0.25 * exp( -0.75 * log((re_c-0.25)*(re_c-0.25)+im_c*im_c));
#endif

#ifdef ZPRIME_PRIME_DIVERGENT_PART
         /* --------------------------------------------------------- */
         /* the interesting one is the z-prime-prime */
			/* This one does zpp/N i.e. the normalized divergent part */
         /* here we use a taylor expansion to extrapolate to tau=0 */
			/* This one extrapolates zpp/norm and thus can show only 
			 * divergent term */

         /* first, we need the derivatives of modulus w.r.t tau */
         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         double mp = (sum_ddrep * sum_ddre + sum_ddimp * sum_ddim) / modulus;
         double mpp  = sum_ddrep * sum_ddrep + sum_ddre * sum_ddrepp;
         mpp += sum_ddimp * sum_ddimp + sum_ddim * sum_ddimpp - mp*mp;
         mpp /= modulus;
         
         /* next, we need derivatives of m/n w.r.t tau */
         mod = modulus / sum_n;
         dmod = (mp - mod * sum_np) / sum_n;
         ddmod = (mpp - 2.0 * dmod * sum_np - mod * sum_npp) / sum_n;

         /* finally the taylor expansion */
         glob [i*sizex +j] = mod - tau* (dmod - 0.5 * tau * ddmod);

         // glob [i*sizex +j] -= 0.25 * exp( -0.75 * log((re_c-0.25)*(re_c-0.25)+im_c*im_c));
#endif


#ifdef PLAIN_OLD_Z_NORMALIZED_DIVERGENCE_FREE
         /* --------------------------------------------------------- */
         /* OK, lets do the taylor expansion for just-plain Z */
         /* here we use a taylor expansion to extrapolate to tau=0 */
			/* this takes sum/normalization then subtracts fitted term. 
			 * The resulting image should be identically zero; there should be
			 * nothing left. */

         /* first, we need the drerivatives of modulus w.r.t tau */
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         double mp = (sum_rep * sum_re + sum_imp * sum_im) / modulus;
         double mpp  = sum_rep * sum_rep + sum_re * sum_repp;
         mpp += sum_imp * sum_imp + sum_im * sum_impp - mp*mp;
         mpp /= modulus;
         
         /* next, we need derivatives of m/n w.r.t tau */
         mod = modulus / sum_n;
         dmod = (mp - mod * sum_np) / sum_n;
         ddmod = (mpp - 2.0 * dmod * sum_np - mod * sum_npp) / sum_n;

         /* finally the taylor expansion for the normalized sum */
         glob [i*sizex +j] = mod - tau* (dmod - 0.5 * tau * ddmod);
// printf ("%9.6g	%9.6g	%9.6g\n", re_c, glob[i*sizex+j], mod);

         /* ok, this part should be the divergent part ... */
         theta = 0.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += 0.5;
         tmp = sqrt(re*re+im*im);
         if (0.5 < tmp) tmp = 0.5;

         glob [i*sizex +j] -= tmp;

#endif

// xx #define ZPP_MODULUS_DIVERGENCE_FREE
#ifdef ZPP_MODULUS_DIVERGENCE_FREE
         /* --------------------------------------------------------- */
         /* The interesting one is the z-prime-prime.
			 * This one subtracts divergence from modulus of zpp 
			 * and goes to tau=0.
          * Here, we subtract the leading divergence 
			 * after computing the modulus, not before. 
			 * This is the one which looks to be a modular form of some kind.
			 */

         modulus = sqrt (sum_ddre*sum_ddre + sum_ddim*sum_ddim);
         double mp = (sum_ddrep * sum_ddre + sum_ddimp * sum_ddim) / modulus;
         double mpp  = sum_ddrep * sum_ddrep + sum_ddre * sum_ddrepp;
         mpp += sum_ddimp * sum_ddimp + sum_ddim * sum_ddimpp - mp*mp;
         mpp /= modulus;
         
         /* finally the taylor expansion */
         /* subtract the main-body divergence */
         tmp = 0.25 * exp( -0.75 * log((re_c-0.25)*(re_c-0.25)+im_c*im_c));
         
         glob [i*sizex +j] = (modulus-tmp*sum_n);
         glob [i*sizex +j] -= tau* ((mp-tmp*sum_np) - 0.5 * tau * (mpp-tmp*sum_npp));
//glob [i*sizex +j] *= variance;
//glob [i*sizex +j] = variance;
#endif

// #define COMPLEX_ZPP_MINUS_DIVERGENCE
#ifdef COMPLEX_ZPP_MINUS_DIVERGENCE
         /* --------------------------------------------------------- */
         /* The interesting one is the z-prime-prime */
			/* This one subtracts divergence from zpp before taking modulus */

         /* The taylor expansion */
			double zre = sum_ddre - tau *(sum_ddrep - 0.5 * tau *sum_ddrepp);
			double zim = sum_ddim - tau *(sum_ddimp - 0.5 * tau *sum_ddimpp);

			/* Divergence term == 0.25/ (0.25-c)^3/2  */
         double thet = -1.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = 0.25 * pow (mod, -0.75);
         re = mod * cos(thet);
         im = mod * sin(thet);

         zre -= re * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));
         zim -= im * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));

         modulus = sqrt (zre*zre + zim*zim);

         glob [i*sizex +j] = zre;

#define MODULAR_FORM_CORRECTIONS
#ifdef MODULAR_FORM_CORRECTIONS
         glob [i*sizex +j] = zre*re_2v - zim*im_2v;
         glob [i*sizex +j] = zre*re_var - zim*im_var;
         glob [i*sizex +j] = zre*re_iv - zim*im_iv;
         glob [i*sizex +j] = zre*re_sca - zim*im_sca;

         double rem = zre*re_sca - zim*im_sca;
         double imm = zre*im_sca + zim*re_sca;
         modulus = sqrt (rem*rem + imm*imm);
         glob [i*sizex +j] = modulus;
#endif /* MODULAR_FORM_CORRECTIONS */
         glob [i*sizex +j] = modulus;
#endif

// #define PLAIN_OLD_MODULUS_Z_MINUS_DIVERGENCE
#ifdef PLAIN_OLD_MODULUS_Z_MINUS_DIVERGENCE
         /* --------------------------------------------------------- */
         /* OK, lets do the taylor expansion for just-plain modulus of Z */
			/* Perform the tau expanstion to extrapolate */
         /* here, we subtract the leading divergence */
         modulus = sqrt (sum_re*sum_re + sum_im*sum_im);
         double mp = (sum_rep * sum_re + sum_imp * sum_im) / modulus;
         double mpp  = sum_rep * sum_rep + sum_re * sum_repp;
         mpp += sum_imp * sum_imp + sum_im * sum_impp - mp*mp;
         mpp /= modulus;
         
         /* finally the taylor expansion */
         glob [i*sizex +j] = modulus - tau* (mp - 0.5 * tau * mpp);

			/* Divergence term == 1/2 - sqrt (1/4-c)  */
         theta = 0.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += 0.5;
         tmp = sqrt(re*re+im*im);

			/* Fix to make it at 1/2 on the large left bulb */
         if (0.5 < tmp) tmp = 0.5;

         glob [i*sizex +j] -= tmp * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));
#endif

// #define COMPLEX_Z_MINUS_DIVERGENCE
#ifdef COMPLEX_Z_MINUS_DIVERGENCE
         /* --------------------------------------------------------- */
         /* OK, lets do the taylor expansion for just-plain Z in full complex glory */
			/* That is do it for z and not for the modulus */
			/* Perform the tau expanstion to extrapolate */
         /* here, we subtract the leading divergence */
         
         /* Now the taylor expansion */
			double zre = sum_re - tau *(sum_rep - 0.5 * tau *sum_repp);
			double zim = sum_im - tau *(sum_imp - 0.5 * tau *sum_impp);

			/* Divergence term == 1/2 - sqrt (1/4-c)  */
         theta = 0.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);
         re += 0.5;

			/* Fix to make it at 1/2 on the large left bulb */
         tmp = sqrt(re*re+im*im);
         if (0.5 < tmp) { re = -0.5; im = 0.0; }

         zre -= re * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));
         zim -= im * (sum_n - tau* (sum_np - 0.5 * tau * sum_npp));

         modulus = sqrt (zre*zre + zim*zim);

         glob [i*sizex +j] = modulus;
#endif

#if WHATEVER
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

         theta = -1.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = pow (mod, -0.75);
         re = 0.25 * mod * cos(theta);
         im = 0.25 * mod * sin(theta);

         glob [i*sizex +j] =  (sum_ddre/sum_n-re)*(sum_ddre/sum_n-re);
         glob [i*sizex +j] += (sum_ddim/sum_n-im)*(sum_ddim/sum_n-im);
         glob [i*sizex +j] = sqrt (glob[i*sizex+j]);

         /* --------------------------------------------------------- */
         theta = 0.5 * atan2 (-im_c, 0.25-re_c);
         mod = (re_c-0.25)*(re_c-0.25)+im_c*im_c;
         mod = pow (mod, 0.25);
         re = - mod * cos(theta);
         im = - mod * sin(theta);

         re += re_c;
         im += im_c;

         glob [i*sizex +j] =  sqrt (re*re +im*im);
#endif
         /* --------------------------------------------------------- */
         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
