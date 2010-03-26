/*
 * brat.C
 *
 * FUNCTION:
 * Explore Hausdorf measure of mandelbrot set.
 * And other stuff.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "brat.h"
#include "Farey.h"
#include "FareyTree.h"

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_out (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax)
{
   int		i,j, globlen, itermax_orig;
   double	re_start, im_start, deltax, deltay;
   double	re_position, im_position;
   double	re_c, im_c;
   double	re, im, tmp;
   int		loop;
   double modulus=0.0, frac, mu;
   double escape_radius = 3.1;
   double ren, tl;
   // double r, phi;

   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0/ log(2.0);
   itermax_orig = itermax;
   

   deltax = width / (double) sizex;
   deltay = height / (double) sizey;
   re_start = re_center - width / 2.0;
   im_start = im_center + height / 2.0;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_c = re_position;
         im_c = im_position;

         /* mapping tricks */
         // re_c = re_position - (re_position*re_position-im_position*im_position);
         // im_c = im_position - 2.0 * re_position * im_position;

#ifdef FLATTEN_CARDIOD_MAP

         /* map to cardiod lam(1-lam) */
         r = im_position;
         phi = M_PI*re_position;

         phi = 2.0*M_PI / re_position;
          
         /*
         works pretty well
         r *= r;
         r = 1.0 + (r-1.0)*sin(0.5*phi)*sin(0.5*phi);
         */
         r -= 1.0;
         // website was done with sin^2(phi/2) but 1-cos(phi/2) is
         // better, according to "Makc the great"
         // r *= sin(0.5*phi)*sin(0.5*phi);
         r *= 1.0 - cos(0.5*phi);
         // r *= (0.5*phi)*sin(0.5*phi);
         // r *= (0.5*phi)* (0.5*phi);
        
         r += 1.0;
         re_c = 0.5 * r * (cos (phi) - 0.5 * r * cos (2.0*phi));
         im_c = 0.5 * r * (sin (phi) - 0.5 * r * sin (2.0*phi));

// hack
// itermax = itermax_orig + 5*re_position*re_position;
#endif /* FLATTEN_CARDIOD_MAP */

         /* remaps */
// #define REMAP_CARDIOD
#ifdef REMAP_CARDIOD
         {
         double u,v;

         re = 0.25 - re_position;
         im = -im_position;
         u = sqrt (0.5*(re + sqrt (re*re + im*im)));
         v = 0.5 * im / u;
         u = 0.5 - u;
         r = sqrt (u*u + v*v);
         phi = atan2 (v,u);
         if (0.0 > phi) phi += 2.0*M_PI;
         phi /= 2.0*M_PI;
         /* phi runs from 0 to 1 */

         phi = phi / (1.0+phi);

         r = 0.5 *sqrt(2.0*r);
         // r = 0.5*(r+0.5);
         
         phi *= 2.0*M_PI;
         re_c = r * (cos (phi) - r * cos (2.0*phi));
         im_c = r * (sin (phi) - r * sin (2.0*phi));

         // printf ("start (%f %f) map to (%f %f)\n", re_position, im_position, re_c, im_c);
         }
#endif  /* REMAP_CARDIOD */

         re = re_c;
         im = im_c;
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;
// printf ("its %d %g %g \n", loop, re, im);
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         modulus = sqrt (modulus);
         frac = log (log (modulus)) *tl;
         mu = ((double) (loop+1)) - frac;

         /* frac =  Re (c/z*z) */
/*
         tmp = (re*re-im*im);
         im = 2.0*re*im;
         re = tmp;
         frac = re*re_c + im*im_c;
         frac /= re*re+im*im;
         if (0.0 > frac) frac = -frac;
*/
         
#ifdef JUNK
         /* second order correction */
         frac = frac * ( 2.0 - frac + 2.0*ren);
       
         nprev = loop - nprev;
         prev = frac - prev - (double)nprev;
/*
         printf ("ij %d %d deltan %d delta frac %f \n", i, j, nprev, prev);
*/
         nprev = loop;
         prev = frac;
         glob [i*sizex +j] = ((double) loop) - frac; 
#endif 

         glob [i*sizex +j] = frac; 
         if (loop<itermax) {
            glob [i*sizex +j] = sqrt(sqrt ((((double) loop) -frac)/ ((double) itermax))); 
         } else {
            glob [i*sizex +j] = 0.0;
         }

/*
         glob [i*sizex +j] = ((float) (loop%10)) / 10.0; 
if (loop == itermax) {
glob[i*sizex+j] = 0.0; } else {glob[i*sizex+j]=0.9999;}
*/

/*
if (0.25*width*width > ((re_position-re_center)* (re_position-re_center) +
(im_position-im_center)* (im_position-im_center))){
         glob [i*sizex +j] -= 0.3;
} else {
         glob [i*sizex +j] += 0.3;
} 
*/

         re_position += deltax;
      }
      im_position -= deltay;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* utility function performs iteration, returns results of iteration */

#define OFL  (1e180)
#define DOFL  (1e20)
#define DMAX  (1e50)

void 
iterate (double re_c,
         double im_c,
         double escape_radius,
         double &re,
         double &im,
         double &dre,
         double &dim,
         int &itermax
         ) 
{
   int loop;
   double esq, tmp, modulus;

   esq = escape_radius*escape_radius;

   re = 0.0;
   im = 0.0;
   dre = 0.0;
   dim = 0.0;
   for (loop=1; loop <=itermax; loop++) {
      /* compute infinitessimal flow */
      tmp = 2.0 * (re*dre - im*dim) +1.0;
      dim = 2.0 * (re*dim + im*dre);
      dre = tmp;

      /* compute actual */
      tmp = re*re - im*im + re_c;
      im = 2.0*re*im + im_c;
      re = tmp;

      while ((OFL < re) || (OFL < im)) {
         re /= DOFL;
         im /= DOFL;
         dre /= DOFL;
         dim /= DOFL;
      }


      modulus = (re*re + im*im);
      if (modulus > esq) {loop++; break; }
   }    

   /* in case we itermaxed, subtract 1 */
   itermax = loop-1;  /* actual count */
}


/*-------------------------------------------------------------------*/
/* regulated mandelbrot exterior, terminates
 * after a finite number of iterations.
 */

void mandelbrot_out_regulated (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
	double regulator)
{
   int		i,j, globlen, itermax_orig;
   double	re_start, im_start, deltax, deltay;
   double	re_position, im_position;
   double	re_c, im_c;
   double	re, im, tmp;
   int		loop;
   double modulus=0.0, frac, mu;
   double escape_radius = 3.1, esq;
   double ren, tl;

   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0/ log(2.0);
   itermax_orig = itermax;
   
	regulator = exp (- regulator);

   deltax = width / (double) sizex;
   deltay = height / (double) sizey;
   re_start = re_center - width / 2.0;
   im_start = im_center + height / 2.0;
	esq = escape_radius * escape_radius;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_c = re_position;
         im_c = im_position;

         re = re_c;
         im = im_c;
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > esq) break;

            im_c *= regulator;
            re_c *= regulator;
         }    

         modulus = sqrt (modulus);
         frac = log (log (modulus)) *tl;
         mu = ((double) (loop+1)) - frac;

         glob [i*sizex +j] = frac; 
         if (loop<itermax) {
            glob [i*sizex +j] = sqrt(sqrt ((((double) loop) -frac)/ ((double) itermax))); 
         } else {
            glob [i*sizex +j] = 0.0;
         }


         re_position += deltax;
      }
      im_position -= deltay;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* This routine fills in the exterior of the mandelbrot set with
 * a shoddy aproximation to the phase.  No good for rays.
 */

void mandelbrot_wind (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax)
{
   int		i,j,k,ii,  il, globlen, itermax_orig;
   double	re_start, im_start, deltax, deltay;
   double	re_position, im_position;
   double	re_c, im_c;
   double	re, im, tmp;
   double 	dre, dim;
   double 	dmu_re, dmu_im;
   int		loop;
   int		*bits;
   double modulus=0.0, frac, mu;
   double ren, otl;
   double phi=0.0, tphi;

   /* the qualty of the results depends on using as large an 
    * escape radius as possible.
    */
   double escape_radius = 1.131e46;

   ren = log( log (escape_radius)) / log(2.0);
   otl = 1.0/ log(2.0);

   itermax --;
   itermax_orig = itermax;

   /* allocate storage for binary tree */
   bits = (int *) malloc ((itermax+50) * sizeof (int));
   

   deltax = width / (double) sizex;
   deltay = height / (double) sizey;
   re_start = re_center - width / 2.0;
   im_start = im_center + height / 2.0;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_c = re_position;
         im_c = im_position;


         /* iterate */
         loop = itermax;
         iterate (re_c, im_c, escape_radius,
                 re, im, dre, dim, loop);
   
         phi = 0.0;
         tphi = 0.0;
         if (loop < itermax) {
            /* compute fractional iteration */
            modulus = (re*re + im*im);
            modulus = sqrt (modulus);
            frac = log (log (modulus)) *otl;
            mu = ((double) (loop+1)) - frac;
   
            /* compute the derivative */
            dmu_re = re*dre + im*dim;
            dmu_im = re*dim - im*dre;
            modulus = (re*re + im*im);
            dmu_re *= -1.0 / modulus;
            dmu_im *= -1.0 / modulus;
            modulus = 1.0 / log(sqrt (modulus));
            dmu_re *= modulus;
            dmu_im *= modulus;
            modulus = 1.0/(dmu_re*dmu_re + dmu_im*dmu_im);
   
            for (k=loop; k>0; k--) 
            {
               /* save up angle */
               /* extract the binary bit */
               phi = atan2(im, re);
               if (0.0<phi) { bits[k] = 0; } else { bits[k] =1; }
               if (0.0>phi) phi+= 2.0*M_PI;
      
               /* redirection for the next step */
               /* uhh, remember its contra not covarient */
// printf ("start %d c=(%g %g) p=%g	z=(%g %g) d=(%g %g)	dm=(%g %g)\n\n", loop,
//  re_c, im_c, phi, re, im, dre, dim, -dmu_re*modulus , dmu_im*modulus);
               re_c -= dmu_re*modulus;
               im_c += dmu_im*modulus;
   
   
               /* lets try fitting */
               tphi = -1000.0;
               for (ii=0; ii<1; ii++)
               { 
                  double gmu, gphi;
                  double dt_re, dt_im;
                  int lp;
      
                  /* ok, now try it */
                  lp = k-1;
                  iterate (re_c, im_c, DMAX*escape_radius,
                       re, im, dre, dim, lp);
      
                  /* compute fractional iteration */
                  modulus = (re*re + im*im);
                  modulus = sqrt (modulus);
                  frac = log (log (modulus)) *otl;
                  gmu = ((double) (lp+1)) +1.0 - frac;
      
                  /* guesstimate angle */
                  gphi = atan2(im, re);
                  if (0.0>gphi) gphi+= 2.0*M_PI;
      
                  if (-100.0>tphi) {
                     tphi = 0.5*phi;
                     if (M_PI<gphi) tphi += M_PI;
                  }
      
                  /* compute the derivative */
                  dmu_re = re*dre + im*dim;
                  dmu_im = re*dim - im*dre;
                  modulus = (re*re + im*im);
                  dmu_re *= 1.0 / modulus;
                  dmu_im *= 1.0 / modulus;
                  dt_re = -dmu_im;
                  dt_im = dmu_re;
                  modulus = - 1.0 / log(sqrt (modulus));
                  dmu_re *= modulus;
                  dmu_im *= modulus;

// printf ("its %d %d c=(%g %g)	gp=%g p=%g  gm=%g m=%g	dm=(%g %g)\n", ii, lp,
//   re_c, im_c, gphi, tphi, gmu, mu, dt_re , dt_im);

                  /* refine the guess */
                  modulus = 1.0/(dmu_re*dmu_re + dmu_im*dmu_im);
#if 0
                  re_c -= (gmu-mu)*dmu_re*modulus;
                  im_c += (gmu-mu)*dmu_im*modulus;
      
                  re_c += (gphi-tphi)*dt_re;
                  im_c += (gphi-tphi)*dt_im;
#endif
               }
// printf ("\n");

   
            }

   
            // binary construction of angle
            tphi = 0.0;
            tmp = 1.0;
            // go no more an 2**48 in the number of bits 
            il = (loop>48) ? 48:loop;
            for (k=1; k<il; k++) {
               tmp *= 0.5;
               if (bits[k]) tphi += tmp;
            }
            // phi = tphi;
   
         }

         phi = 0.5*phi/M_PI;

         glob [i*sizex +j] = phi;

         re_position += deltax;
      }
      im_position -= deltay;  /*top to bottom, not bottom to top */
   }

   free (bits);
}

/*-------------------------------------------------------------------*/
/* This routine fills in the exterior of the mandelbrot set with
 * the phase/ winding number.  Has a few discontinuities/imperfections,
 * but is othrewise excellent at drawing rays.  i.e. phase is exactly 
 * correct.
 */

void mandelbrot_windsimple (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   int		nrays,
   double	ppp,
   double 	qqq)
{
   int		i,j,k, globlen, itermax_orig;
   double	re_start, im_start, deltax, deltay;
   double	re_position, im_position;
   double	re_c, im_c;
   double	re, im, tmp;
	double   re_prev, im_prev;
	double   re_preprev, im_preprev;
   int		loop;
   double modulus=0.0;
   double escape_radius = 1131.1;
   double ren, otl;
   double phi=0.0, phi_last, phi_c, h_phi_c;
   int wind =0;

   ren = log( log (escape_radius)) / log(2.0);
   otl = 1.0/ log(2.0);

   itermax --;
   itermax_orig = itermax;

   deltax = width / (double) sizex;
   deltay = height / (double) sizey;
   re_start = re_center - width / 2.0;
   im_start = im_center + height / 2.0;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_c = re_position;
         im_c = im_position;

         phi_c = atan2 (im_c, re_c);
         if (0.0 > phi_c) phi_c += 2.0*M_PI;
         h_phi_c = 0.5*phi_c;

         re = 0.0;
         im = 0.0;
         re_prev = 0.0;
			im_prev = 0.0;
         phi_last = -0.01;
         wind = 0;
// printf ("\n phi_c= %12.8g   c=(%g %g)\n", phi_c, re_c, im_c);
         for (loop=1; loop <itermax; loop++) {
				re_preprev = re_prev;
				im_preprev = im_prev;
				re_prev = re;
				im_prev = im;
            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;

            phi = atan2 (im, re);
            // phi = atan2 (im, re+0.5);
            // phi = atan2 (im, re-0.5);
            // phi = atan2 (im-im_prev, re-re_prev);
            // phi = atan2 (im-im_preprev, re-re_preprev);
            // phi = atan2 (2.0*im-im_prev, 2.0*re-re_prev);
            if (0.0>phi) phi += 2.0*M_PI;
            wind += wind;
				
#if WORKS_WITH_SMALL_DEFECTS

				// The below cancels out horns 
            if (phi < phi_last) wind +=2;

// printf ("n=%d %12.8g %12.8g %12.8g %12.8g %10.6g %d\n", 
// loop, phi, phi_last, 2.0*phi_last-phi, frac,
// phi+((double)wind)*M_PI, wind);

				/* The following is critical for going from discontinuous
				 * power-two bars to smooth func. */
            /* if northern half else southern half */
            if (M_PI > phi_c) { 
               while (M_PI < phi) {
                  phi -= M_PI;
                  wind ++;
               }
            } else {
               while (M_PI > phi) {
                  phi += M_PI;
                  wind --;
               }
            }
#endif /* WORKS_WITH_SMALL_DEFECTS */

#define THEORY
#ifdef THEORY
				phi_last *= 2.0;
				if (2.0 * M_PI < phi_last) wind ++;

				if ((M_PI < phi_last) && (2.0 * M_PI > phi_last))
				{
					if ((0 < phi) && (M_PI > phi)) wind ++;
				}
#endif /* THEORY */
				

            phi_last = phi;

            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         phi /= 2.0*M_PI;

         // use the winding number count to compute the smooth phase
         if (loop < itermax) {
            phi += 0.5 * ((double) wind);

            // the phase winds around 2pi *2^(loop-1) times
            phi /= pow (2.0, (double) (loop-1));

         }

         // colorize the landing rays
         if (0 < nrays) {
            double tphi = phi * 2.0 * ((double) nrays);
   
            k = (int) tphi;
            if (k%2) { 
               tphi -= (double)k;
            } else {
               tphi = (double)(k+1) -tphi;
            }

            if (phi < ppp/((double)nrays)) tphi *= 0.7;
            if (phi > qqq/((double)nrays)) tphi *= 0.7;

            tphi *= tphi;
            tphi *= tphi;
            phi = tphi;
         }

         glob [i*sizex +j] = phi;

         re_position += deltax;
      }
      im_position -= deltay;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_inverse (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   int ire, jim;
   double modulus=0.0, phase, tp, frac;
   double escape_radius = 2000.1;
   double ren, tl;
   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0/ log(2.0);
   

   delta = width / (double) sizex;
delta *= 0.25;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<4*sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<4*sizex; j++) {
         re = re_position;
         im = im_position;
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         if (loop == itermax) continue;

         modulus = sqrt (modulus);
         frac = ren - log (log (modulus)) *tl;
         frac = ((double) loop) - ren + frac; 

         modulus = log (modulus);
         phase = atan2 (im, re);
         tp = 1.0 / pow (2.0, loop);
         modulus *= tp;
         phase *= tp;
         modulus = exp (modulus);
         re = cos (phase) * modulus;
         im = sin (phase) * modulus;
         
         re -= 1.0;

         im /= re;

         re *= 24.0;
         im *= 1.0;

         ire = (int) (sizex * (re + 0.5));
         if (0 > ire) ire = 0;
         if (sizex <= ire) ire = sizex-1;

         im *= ((double) sizey)/ ((double) sizex);
         jim = (int) (sizey * (im + 0.5));

         if (0 > jim) jim = 0;
         if (sizey <= jim) jim = sizey-1;

         glob [jim*sizex +ire] ++;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_stop (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re = re_position;
         im = im_position;
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            if ((re*re + im*im) > 154.0) break;
         }    
         /* glob [i*sizex +j] = (2.0+ re)/2.5;  */
         /* glob [i*sizex +j] = re; */
         glob [i*sizex +j] = im;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the Cliff Pickover's stalks algorithm */

void mandelbrot_stalk (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double	stalkx,
   double	stalky)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re, im, tmp, tmpx, tmpy;
   // double	visited_x, visited_y;
   int		loop;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re = re_position;
         im = im_position;
#ifdef MAXDIST
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            if ((re*re + im*im) > 154.0) break;
         }
         visited_x = re;
         visited_y = im;
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
#endif // MAXDIST
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            if ((re*re + im*im) > 154.0) break;

#ifdef HAND_BUILT_COLORMAP
            Hand-built colormap -- visually bad idea
            if (((-0.5 < re) && (re < 0.5)) ||((-0.5 < im) && (im < 0.5))) 
                       if (0.1 > glob [i*sizex +j]) glob [i*sizex +j] = 0.1;
            if (((-0.2 < re) && (re < 0.2)) ||((-0.2 < im) && (im < 0.2))) 
                       if (0.2 > glob [i*sizex +j]) glob [i*sizex +j] = 0.2;
            if (((-0.1 < re) && (re < 0.1)) ||((-0.1 < im) && (im < 0.1))) 
                       if (0.3 > glob [i*sizex +j]) glob [i*sizex +j] = 0.3;
            if (((-0.05 < re) && (re < 0.05)) ||((-0.05 < im) && (im < 0.05))) 
                       if (0.5 > glob [i*sizex +j]) glob [i*sizex +j] = 0.5;
            if (((-0.02 < re) && (re < 0.02)) ||((-0.02 < im) && (im < 0.02))) 
                       if (0.6 > glob [i*sizex +j]) glob [i*sizex +j] = 0.6;
            if (((-0.01 < re) && (re < 0.01)) ||((-0.01 < im) && (im < 0.01))) 
                       if (0.7 > glob [i*sizex +j]) glob [i*sizex +j] = 0.7;
            if (((-0.005 < re) && (re < 0.005)) ||((-0.005 < im) && (im < 0.005))) 
                       if (0.8 > glob [i*sizex +j]) glob [i*sizex +j] = 0.8;
            if (((-0.002 < re) && (re < 0.002)) ||((-0.002 < im) && (im < 0.002))) 
                       if (0.9 > glob [i*sizex +j]) glob [i*sizex +j] = 0.9;
            if (((-0.001 < re) && (re < 0.001)) ||((-0.001 < im) && (im < 0.001))) 
                       if (1.0 > glob [i*sizex +j]) glob [i*sizex +j] = 1.0;
#endif /* HAND_BUILT_COLORMAP */

#ifdef STALK
            tmpx = re-stalkx;
            if (0.0 > tmpx) tmpx = -tmpx;
            tmpy = im-stalky;
            if (0.0 > tmpy) tmpy = -tmpy;

            ltmpx = -log (tmpx); 
            if (ltmpx > glob [i*sizex +j]) glob [i*sizex +j] = ltmpx;

            ltmpy = -log (tmpy); 
            if (ltmpy > glob [i*sizex +j]) glob [i*sizex +j] = ltmpy;
#endif /* STALK */

#ifdef CROSS
            tmpx = re-stalkx;
            if (0.0 > tmpx) tmpx = -tmpx;
            tmpy = im-stalky;
            if (0.0 > tmpy) tmpy = -tmpy;

            ltmpx = -log (tmpx); 
            if (0.2 > tmpy) {
               if (ltmpx > glob [i*sizex +j]) glob [i*sizex +j] = ltmpx;
            }

            ltmpy = -log (tmpy); 
            if (0.2 > tmpx) {
               if (ltmpy > glob [i*sizex +j]) glob [i*sizex +j] = ltmpy;
            }
#endif /* CROSS */

#define CIRCSTALK
#ifdef CIRCSTALK
            tmpx = re - stalkx;
            tmpy = im - stalky;
            tmp = tmpx*tmpx + tmpy*tmpy -0.04;
            tmp = tmpx*tmpx + tmpy*tmpy;
            if (0.0 > tmp) tmp = -tmp;

            tmp = -log (tmp); 
            if (tmp > glob [i*sizex +j]) glob [i*sizex +j] = tmp;
#endif // CIRCSTALK

           
#ifdef MAXDIST
            tmpx = re-stalkx - visited_x;
            tmpy = im-stalky - visited_y;

            tmpx = re - stalkx;
            tmpy = im - stalky;
            tmp = tmpx*tmpx + tmpy*tmpy;

            tmp = -log (tmp); 
            tmpx = re - stalkx-0.001;
            tmpy = im - stalky-0.001;

            tmp += log (tmpx*tmpx + tmpy*tmpy);

            if (tmp > glob [i*sizex +j]) glob [i*sizex +j] = tmp;
#endif // MAXDIST

           


         }    

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the circle map set using */
/* the classic algorithm */

void circle_out (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	time)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re_omega, im_omega;
   double 	re_K, im_K;
   double	re, im, tmp;
   int		n, loop;
   double	ep, em, es, ec;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.000000001;

/* 
   re_omega = time;
   im_omega = 0.1;
*/

   im_K = time;
   im_omega = 0.0;

   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (0 == i%30) printf(" start row %d of %d\n", i, sizey);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_K = re_position;
         /* im_K = im_position; */
         re_omega = im_position; 
         re = 0.0;
         im = 0.0;

         for (loop=1; loop <itermax; loop++) {
            tmp = 2.0 * M_PI * im;
            ep = exp (tmp);
            em = exp (-tmp);
            es = 0.5 * (ep+em);
            ec = 0.5 * (ep-em);
            tmp = 2.0 * M_PI * re;
            es *= sin (tmp);
            ec *= cos (tmp);

            tmp = re + re_omega - re_K * es + im_K * ec;
            im = im + im_omega - im_K * es - re_K * ec;
            re = tmp;


            if (re*re > 1.0e6) break;
            n = (int) re;
            if (0.0 > re) n--;
            re -= (double) n;

            if (im*im > 2284.0) break;
         }    

         glob [i*sizex +j] = loop; 

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }

   tmp = 1.0 / ((double) itermax);
   for (i=0; i<sizex*sizey; i++) {
      glob[i] *= tmp;
   }
   
}

/*-------------------------------------------------------------------*/
extern double drand48 ();
#define LOOP_COUNT 348 
/* #define LOOP_COUNT 4480 */
#define SETTLE_COUNT 250

#define CBOX_IM_SLOPE 3.0
#define CBOX_IM_CEPT -1.5
#define CBOX_RE_SLOPE 3.0
#define CBOX_RE_CEPT -1.5

/* this routine fills in the interior of the circle map using */
/* the classic algorithm */

void circle_in (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	time)
{
   int	i, globlen;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re_omega, im_omega;
   double 	re_K, im_K;
   double	re, im, tmp;
   int		n, loop;
   double	ep, em, es, ec;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		isamp, irow, horiz_pix, vert_pix;
   double	xs, ys;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.000000001;

   re_omega = time;
   im_omega = 0.1;
   im_omega = 0.0;

/*
   im_K = time;
   im_omega = 0.0;
*/
   isamp = (int) (CBOX_IM_SLOPE*CBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = CBOX_IM_SLOPE * xs + CBOX_IM_CEPT;
      re_position = CBOX_RE_SLOPE * ys + CBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %d of %d\n", i/irow, sizey);
   
      re_K = re_position;
      im_K = im_position; 
      /* re_omega = im_position; */
      re = 0.0;
      im = 0.0;

      for (loop=1; loop <LOOP_COUNT; loop++) {
         tmp = 2.0 * M_PI * im;
         ep = exp (tmp);
         em = exp (-tmp);
         es = 0.5 * (ep+em);
         ec = 0.5 * (ep-em);
         tmp = 2.0 * M_PI * re;
         es *= sin (tmp);
         ec *= cos (tmp);

         tmp = re + re_omega - re_K * es + im_K * ec;
         im = im + im_omega - im_K * es - re_K * ec;
         re = tmp;


         if (re*re > 1.0e6) break;
         n = (int) re;
         if (0.0 > re) n--;
         re -= (double) n;

         if (im*im > 2284.0) break;

         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
             (vert_pix >= 0) && (vert_pix < sizey) &&
             (SETTLE_COUNT < loop)) {
             glob [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   tmp = 1.0 / ((double) itermax*(LOOP_COUNT-SETTLE_COUNT)) ;
   for (i=0; i<sizex*sizey; i++) {
      glob[i] *= tmp;
   }
   
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */

#define BBOX_IM_SLOPE 3.0
#define BBOX_IM_CEPT -1.5
#define BBOX_RE_SLOPE 3.0
#define BBOX_RE_CEPT -2.1

void 
mandelbrot_measure (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.00001;
   
   /* random seeds start (ideally) with uniform density inside the 
    * mandelbrot set.  With a slight loss of efficiency, and no loss of
    * correctness, it is sufficient to sart with a uniform density of 
    * points within the bounding box of the mandelbrot set.
    *
    * Since most of the random hits will occur outside of the region
    * we are interested in displaying, we need to increase the number
    * of random starts in inverse proportion of the size of the window
    * to the size of the mandelbrot bounding box.
    *
    * To get any reasonable sort of sampling, we need the number of 
    * random seeds to be proportional to the number of pixels.
    * 
    * To make sure we get into some reasonably stable orbits,
    * must iterate each seed at least 20 times (this is an untested
    * hypothesis)
    */
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < SETTLE_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
      }    

      for (loop=0; loop < LOOP_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
            glob [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {
      glob [i] = glob[i] / (double) (renorm*itermax*LOOP_COUNT) ;
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */

void 
mandelbrot_offset (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 1.0e-10;
   
   /* random seeds start (ideally) with uniform density inside the 
    * mandelbrot set.  With a slight loss of efficiency, and no loss of
    * correctness, it is sufficient to sart with a uniform density of 
    * points within the bounding box of the mandelbrot set.
    *
    * Since most of the random hits will occur outside of the region
    * we are interested in displaying, we need to increase the number
    * of random starts in inverse proportion of the size of the window
    * to the size of the mandelbrot bounding box.
    *
    * To get any reasonable sort of sampling, we need the number of 
    * random seeds to be proportional to the number of pixels.
    * 
    * To make sure we get into some reasonably stable orbits,
    * must iterate each seed at least 20 times (this is an untested
    * hypothesis)
    */
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < SETTLE_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
      }    

      for (loop=0; loop < LOOP_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * (re-re_position) - x_off);
         vert_pix = (int) (y_slope * (im-im_position) - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
            glob [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {
      glob [i] = glob[i] / (double) (itermax*LOOP_COUNT) ;
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. It keeps track of two things: the average age, 
 * and the mean square age deviation. */

void 
mandelbrot_age (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float	*count;
   float	*square, *cube, *quart, *quint;
   double ollie;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   square = (float *)  malloc (globlen * sizeof (float));
   cube = (float *)  malloc (globlen * sizeof (float));
   quart = (float *)  malloc (globlen * sizeof (float));
   quint = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 1.0e-10;
      count [i] = 1.0e-10;
      square [i] = 1.0e-10;
      cube [i] = 1.0e-10;
      quart [i] = 1.0e-10;
      quint [i] = 1.0e-10;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
            glob [vert_pix*sizex + horiz_pix] += loop;
            square [vert_pix*sizex + horiz_pix] += loop*loop;
            cube [vert_pix*sizex + horiz_pix] += loop*loop*loop;
            quart [vert_pix*sizex + horiz_pix] += loop*loop*loop*loop;
            quint [vert_pix*sizex + horiz_pix] += loop*loop*loop*loop*loop;
            count [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   ollie = 1.0 / ((double) LOOP_COUNT);

   for (i=0; i<globlen; i++) {
      double tmp;

#ifdef FIRST_MOM
      /* compute average age */
      tmp = 1.0 / count[i];
      glob [i] = glob[i] * tmp * ollie;
#endif 
    
#ifdef SECOND_MOM
      /* compute mean square age */
      tmp = 1.0 / count[i];
      glob[i] = tmp *(square[i] - glob[i]*glob[i] *tmp);
      glob[i] *= ollie*ollie;
#endif

#ifdef THIRD_MOM
      /* compute third moment */
      tmp = 1.0 / count[i];
      glob[i] *= tmp;
      square[i] *= tmp*tmp;
      cube[i] *= tmp*tmp*tmp;
      glob[i] = cube[i] - 3.0*square[i]*glob[i] + 2.0*glob[i]*glob[i]*glob[i];
      glob[i] *= ollie*ollie*ollie;
#endif

#ifdef FOURTH_MOM
      /* compute third moment */
      tmp = 1.0 / count[i];
      glob[i] *= tmp;
      square[i] *= tmp*tmp;
      cube[i] *= tmp*tmp*tmp;
      quart[i] *= tmp*tmp*tmp*tmp;

      tmp = 6.0 * square[i] - 3.0 * glob[i]*glob[i];
      tmp = -4.0 * cube[i] + tmp*glob[i];
      tmp = quart[i] + tmp*glob[i];
      glob[i] = tmp *ollie*ollie*ollie*ollie;
#endif

#define FIFTH_MOM 
#ifdef FIFTH_MOM
      /* compute third moment */
      tmp = 1.0 / count[i];
      glob[i] *= tmp;
      square[i] *= tmp*tmp;
      cube[i] *= tmp*tmp*tmp;
      quart[i] *= tmp*tmp*tmp*tmp;
      quart[i] *= tmp*tmp*tmp*tmp*tmp;

      tmp = -10.0 * square[i] + 4.0 * glob[i]*glob[i];
      tmp = 10.0 * cube[i] + tmp*glob[i];
      tmp = -5.0 * quart[i] + tmp*glob[i];
      tmp = quint[i] + tmp*glob[i];
      glob[i] = tmp *ollie*ollie*ollie*ollie*ollie;
#endif

   }

   free (count);
   free (square);
   free (cube);
   free (quart);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. It keeps track of a rough lyapunov exponenet */

void 
mandelbrot_lyapunov (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i, j;
   long		isamp, irow;
   double	re_start, im_start, delta;
   double	re_position[5], im_position[5];
   double	re[5], im[5], tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float	*count;
   double 	dist, logdelta;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;

   delta = 0.1 * width / ((double) sizex);
   logdelta = log (2.0 * delta);
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 1e-20;
      count [i] = 1e-20;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position[0] = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position[0] = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;

      im_position[1] = im_position[0] + delta;
      re_position[1] = re_position[0];
   
      im_position[2] = im_position[0];
      re_position[2] = re_position[0] + delta;
   
      im_position[3] = im_position[0] - delta;
      re_position[3] = re_position[0];
   
      im_position[4] = im_position[0];
      re_position[4] = re_position[0] - delta;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      for (j=0; j<5; j++) {
         re[j] = im[j] = 0.0;
      }
      for (loop=1; loop < LOOP_COUNT; loop++) {
         for (j=0; j<5; j++) {
            tmp = re[j]*re[j] - im[j]*im[j] + re_position[j];
            im[j] = 2.0*re[j]*im[j] + im_position[j];
            re[j] = tmp;
         }
         if ((re[0]*re[0] + im[0]*im[0]) > 7.0) break;
         if ((re[1]*re[1] + im[1]*im[1]) > 7.0) break;
         if ((re[2]*re[2] + im[2]*im[2]) > 7.0) break;
         if ((re[3]*re[3] + im[3]*im[3]) > 7.0) break;
         if ((re[4]*re[4] + im[4]*im[4]) > 7.0) break;
#ifdef STRAIGHT
         horiz_pix = (int) (x_slope * re[0] - x_off);
         vert_pix = (int) (y_slope * im[0] - y_off);
#else
         horiz_pix = (int) (x_slope * (re[0]-re_position[0]) - x_off);
         vert_pix = (int) (y_slope * (im[0]-im_position[0]) - y_off);
#endif
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {

            dist = 0.0;
            for (j=1; j<5; j++) {
               tmp = im[j] - im[0]; 
               dist += tmp*tmp;
               tmp = re[j] - re[0]; 
               dist += tmp*tmp;
            }
            dist = sqrt (dist);
            dist = log (dist) - logdelta;
            dist /= (double) loop;

            glob [vert_pix*sizex + horiz_pix] += dist;
            count [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {

      /* compute average age */
      glob [i] = glob[i] / (double) (LOOP_COUNT*count[i]) ;
    
   }

   free (count);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. It keeps track of a migration direction */

void 
mandelbrot_migrate (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i, j;
   long		isamp, irow;
   double	re_start, im_start, delta;
   double	re_position[5], im_position[5];
   double	re[5], im[5], tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float	*count;
   double 	logdelta;
   float 	*vecx, *vecy;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;

   delta = 0.1 * width / ((double) sizex);
   logdelta = log (2.0 * delta);
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   vecx = (float *)  malloc (globlen * sizeof (float));
   vecy = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 1e-20;
      count [i] = 1e-20;
      vecx [i] = 1e-20;
      vecy [i] = 1e-20;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position[0] = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position[0] = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;

      im_position[1] = im_position[0] + delta;
      re_position[1] = re_position[0];
   
      im_position[2] = im_position[0];
      re_position[2] = re_position[0] + delta;
   
      im_position[3] = im_position[0] - delta;
      re_position[3] = re_position[0];
   
      im_position[4] = im_position[0];
      re_position[4] = re_position[0] - delta;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      for (j=0; j<5; j++) {
         re[j] = im[j] = 0.0;
      }
      for (loop=1; loop < LOOP_COUNT; loop++) {
         for (j=0; j<5; j++) {
            tmp = re[j]*re[j] - im[j]*im[j] + re_position[j];
            im[j] = 2.0*re[j]*im[j] + im_position[j];
            re[j] = tmp;
         }
         if ((re[0]*re[0] + im[0]*im[0]) > 7.0) break;
         if ((re[1]*re[1] + im[1]*im[1]) > 7.0) break;
         if ((re[2]*re[2] + im[2]*im[2]) > 7.0) break;
         if ((re[3]*re[3] + im[3]*im[3]) > 7.0) break;
         if ((re[4]*re[4] + im[4]*im[4]) > 7.0) break;
         horiz_pix = (int) (x_slope * re[0] - x_off);
         vert_pix = (int) (y_slope * im[0] - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {

            for (j=1; j<5; j++) {
               tmp = (re[j] - re[0]) * (re[j] - re[0]); 
               tmp += (im[j] - im[0]) * (im[j] - im[0]); 
               tmp = 1.0 / sqrt (tmp);
               vecx [vert_pix*sizex + horiz_pix] += (re[j] - re[0]) *tmp; 
               vecy [vert_pix*sizex + horiz_pix] += (im[j] - im[0]) *tmp; 
            }

            count [vert_pix*sizex + horiz_pix] ++;
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {

      /* compute average angle */
      glob [i] = atan2 (vecy [i], vecx[i]); 
      glob [i] = (glob[i] +M_PI)/ (2.0 * M_PI);
   }

   free (vecx);
   free (vecy);
   free (count);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. It keeps track of two things: the average age, 
 * and the mean square age deviation. */

void 
mandelbrot_squige (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float	*count;
   int last;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 0.00001;
      count [i] = 0.00001;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      last = 0;
      re = im = 0.0;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
            glob [vert_pix*sizex + horiz_pix] += last;
            count [vert_pix*sizex + horiz_pix] += 1.0;
            last = (int) count [vert_pix*sizex + horiz_pix];
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {

      /* compute average age */
      glob [i] = glob[i] / (double) (LOOP_COUNT*count[i]) ;
    
   }

   free (count);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */

void 
mandelbrot_phase (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, renew, imnew, dist;
   double	reold, imold;
   double	norm;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float 	*count;
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 0.00001;
      count [i] = 0.00001;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = re_position;
      im = im_position;
      reold = re;
      imold = im;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         renew = re*re - im*im + re_position;
         imnew = 2.0*re*im + im_position;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey) 
             && (2 < loop)) {
/*
            dist = (renew-re)*(renew-re);
            dist += (imnew-im)*(imnew-im);
            dist = atan2 (imnew-im, renew-re);
*/
            dist = (renew-re)*(re-reold) + (imnew-im)*(im-imold);
            norm = (renew-re)*(renew-re) + (imnew-im)*(imnew-im);
            norm *= (re-reold)*(re-reold) + (im-imold)*(im-imold);
            if (1.0e-30 < norm) {
               dist = dist / sqrt (norm);
            } else {
                dist = 0.0;
            }

            glob [vert_pix*sizex + horiz_pix] += dist;
            count [vert_pix*sizex + horiz_pix] ++;
         }

         reold = re;
         imold = im;
         re = renew;
         im = imnew;
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {
      glob [i] = glob[i] / (count[i] * 2.0)  ;
      glob [i] += 0.5;
   }

   free (count);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. The colorations are assigned from the constant
 * term. (the "origin") */

void 
mandelbrot_orig (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, renew, imnew;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float 	*count;
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 0.00001;
      count [i] = 0.00001;
   }
   
   isamp = (int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         renew = re*re - im*im + re_position;
         imnew = 2.0*re*im + im_position;
         if ((re*re + im*im) > 7.0) break;
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
/*
            glob [vert_pix*sizex + horiz_pix] += ys;
*/
            glob [vert_pix*sizex + horiz_pix] = xs;
            count [vert_pix*sizex + horiz_pix] ++;
         }

         re = renew;
         im = imnew;
      }    
   }

   /* renormalize */
/*
   for (i=0; i<globlen; i++) {
      glob [i] = glob[i] / count[i];
   }

*/
   free (count);
}

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm. The colorations are assigned from the 
 * the most recent location of the pixel */

void 
mandelbrot_next (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, renew, imnew;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float 	*count;
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 0.00001;
      count [i] = 0.00001;
   }
   
   isamp = (long int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         renew = re*re - im*im + re_position;
         imnew = 2.0*re*im + im_position;
         if ((re*re + im*im) > 7.0) break;

#define OM 10.0
         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
             xs = sin (OM * renew);
             ys = sin (OM * imnew);

             glob [vert_pix*sizex + horiz_pix] += xs*xs*ys*ys;
             count [vert_pix*sizex + horiz_pix] ++;
             
         }    

         re = renew;
         im = imnew;

      }
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {
      glob [i] = glob[i] / count[i];
   }

   free (count);
}

/*-------------------------------------------------------------------*/

void 
whack (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	renorm)
{
   int		globlen;
   long		i;
   long		isamp, irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   double	xs, ys;
   float	*count;
   float 	ico;
   int		hopix, vopix, loco;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   

   globlen = sizex*sizey;
   count = (float *)  malloc (globlen * sizeof (float));
   for (i=0; i<globlen; i++) {
      glob [i] = 0.00001;
      count [i] = 0.00001;
   }
   
   isamp = (long int) (BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width));
   if (0>= isamp) isamp = 1;
   isamp *=  globlen*itermax;
   irow = isamp / sizey;
   
   for (i=0; i<isamp; i++) {
      xs = drand48();
      ys = drand48();
      im_position = BBOX_IM_SLOPE * xs + BBOX_IM_CEPT;
      re_position = BBOX_RE_SLOPE * ys + BBOX_RE_CEPT;
   
      hopix = (int) (x_slope * re_position - x_off);
      vopix = (int) (y_slope * im_position - y_off);
      if (hopix < 0) hopix = 0;
      if (vopix < 0) vopix = 0;
      if (hopix >= sizex) hopix = sizex-1;
      if (vopix >= sizey) vopix = sizey-1;
      loco = vopix*sizex + hopix;

      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = re_position;
      im = im_position;
      for (loop=0; loop < LOOP_COUNT; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;

         horiz_pix = (int) (x_slope * re - x_off);
         vert_pix = (int) (y_slope * im - y_off);
         if ( (horiz_pix >= 0) && (horiz_pix < sizex) && 
            (vert_pix >= 0) && (vert_pix < sizey)) {
            ico = count [vert_pix*sizex + horiz_pix] ++;

            glob [loco] += ico;
         }
      }    
   }

   /* renormalize */
   for (i=0; i<globlen; i++) {
      /* glob [i] = glob[i] / (double) (LOOP_COUNT*count[i]) ; */
      glob [i] = count[i] / (double) (renorm*itermax*LOOP_COUNT*LOOP_COUNT) ;
   }

   free (count);
}

/*-------------------------------------------------------------------*/
/* This routine fills in the exterior of the mandelbrot set using 
 * an algorithm which is, well, a bit differnet. 
 * I suspect a problem with the original algorithm is that the seed
 * pixels always start on grid boundaries.  The algorithm below
 * provides a random seed instead, and assigns that to the pixel.
 * Like the traditional exterior coloring algorithms, the number of
 * iterations a pixel goes through before it escapes is counted.
 */

void 
random_out (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   int		globlen;
   long		i;
   long		irow;
   double	re_start, im_start;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   double	x_width, y_width;
   double	x_slope, y_slope, x_off, y_off;
   int		horiz_pix, vert_pix;
   
   
   x_width = width;
   y_width = width * ((double) sizey) / ((double) sizex);
   re_start = re_center - x_width / 2.0;
   im_start = im_center - y_width / 2.0;
   x_slope = ((double) sizex) / x_width;
   y_slope = ((double) sizey) / y_width;
   x_off = x_slope * re_start;
   y_off = y_slope * im_start;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   irow = 2*sizey;
   
   for (i=0; i<2*globlen; i++) {
      im_position = y_width * drand48() + im_start;
      re_position = x_width * drand48() + re_start;
   
      if (i%(irow)==0) printf(" start row %ld\n", i/irow);
   
      re = im = 0.0;
      for (loop=0; loop < itermax; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 7.0) break;
      }    

      horiz_pix = (int) (x_slope * re_position - x_off);
      vert_pix = (int) (y_slope * im_position - y_off);
      glob [vert_pix*sizex + horiz_pix] = ((float) (loop%10)) / 10.0; 
   }
}

/*-------------------------------------------------------------------*/
/** This routine does an ordinary height map: the callback is given
 *  an x,y coordinate pair, and is expected to return a single value,
 *  which will be plotted as such.
 */

void 
MakeHeightWrap (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm,
	MakeHeightCB cb)
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
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = cb (re_position, im_position, itermax, renorm);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/** This routine does bifurcation diagrams. The callback is passed a
 * row of pixels, and it is expected to fill out the row, which is then
 * plotted. */

void 
MakeBifurWrap (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm,
	MakeBifurCB cb)
{
   int		i, globlen;
   double	im_start, delta;
   double	im_position;
   
   delta = width / (double) sizex;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);

		cb (&glob[i*sizex], sizex, re_center, width, im_position, itermax, renorm);
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/

extern "C" {
   extern FILE *Fopen(char *name, char *ext);
};

int 
main (int argc, char *argv[]) 
{
   float	*data;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   double	re_center, im_center, width, height;
   int		itermax;
   double	renorm, tmp;
   FILE		*fp;
   int		i, nray=-1;
   double	p=0.0, q=1.0;
	double   funky = 0.0;
   char buff [80];
	char * progname;
   
   if (5 > argc) {
      fprintf (stderr, "Usage: %s <filename> <width> <height> <niter> [<centerx> <centery> <width> [<param>]]\n", argv[0]);
      exit (1);
   }

   itermax = 1;
   data_width = 200;
   data_height = 200;

   data_width = atoi (argv[2]);
   data_height = atoi (argv[3]);
   itermax = atoi (argv[4]);

   data = (float *) malloc (data_width*data_height*sizeof (float));

   re_center = -0.6;
   im_center = 0.0;
   width = 2.8;

   if (6 == argc) {
      funky = atof (argv[5]);
	}
	
   if (8 <= argc) {
      re_center = atof (argv[5]);
      im_center = atof (argv[6]);
      width = atof (argv[7]);
   }
   height = width * ((double) data_height) / ((double) data_width);

   renorm = 1.0;
   if (9 == argc) renorm = atof (argv[8]);

   if (9 < argc) {
      nray = atoi (argv[8]);
      p = atof (argv[9]);
      q = atof (argv[10]);
   }

   printf ("file=%s (%d %d) iter=%d (%f %f) w=%f h=%f param=%f\n", 
        argv[1], data_width, data_height, itermax, 
        re_center, im_center, width, height, renorm); 

   /* Do the interior now */

	progname = argv[0];
	progname = strrchr (progname, '/');
	if (!progname) progname = argv[0];
	else progname ++;

	MakeHisto (data, data_width, data_height,
                  re_center, im_center, width, height, itermax, renorm); 

   if (!strcmp(progname, "brat"))
   mandelbrot_out (data, data_width, data_height,
                  re_center, im_center, width, height, itermax); 
   
   if (!strcmp(progname, "morph"))
   mandelbrot_out_regulated (data, data_width, data_height,
                  re_center, im_center, width, height, itermax, funky); 
   
   if (!strcmp(progname, "wind"))
   mandelbrot_wind (data, data_width, data_height,
                  re_center, im_center, width, height, itermax); 
   
   if (!strcmp(progname, "winds"))
   mandelbrot_windsimple (data, data_width, data_height,
                  re_center, im_center, width, height, itermax, nray, p, q); 
   
   if (!strcmp(progname, "manvert"))
   mandelbrot_inverse (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(progname, "measure"))
   mandelbrot_measure (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(progname, "offset"))
   mandelbrot_offset (data, data_width, data_height,
                  0.5, 0.0, 3.0, itermax, renorm); 
   
   if (!strcmp(progname, "mstop"))
   mandelbrot_stop (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(progname, "stalk"))
   mandelbrot_stalk (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.0, 0.0); 
   
   if (!strcmp(progname, "orig"))
   mandelbrot_orig (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(progname, "phase"))
   mandelbrot_phase (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(progname, "age"))
   mandelbrot_age (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
#ifdef STRAIGHT
   if (!strcmp(progname, "lyapunov"))
   mandelbrot_lyapunov (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
#else
   if (!strcmp(progname, "lyapunov"))
   mandelbrot_lyapunov (data, data_width, data_height,
                  0.5, 0.0, 3.0, itermax, renorm); 
#endif
   

   if (!strcmp(progname, "migrate"))
   mandelbrot_migrate (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(progname, "squige"))
   mandelbrot_squige (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(progname, "next"))
   mandelbrot_next (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(progname, "whack"))
   whack (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   
   if (!strcmp(progname, "circout")) {
      re_center = 0.0;
      im_center = 0.0;
      width = 1.0;
      
      circle_out (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.33); 
   }
   
   if (!strcmp(progname, "circin")) {
      re_center = 0.5;
      im_center = 0.0;
      width = 1.0;
      
      circle_in (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.05); 
   }
   
   /* ---------------------------------------------------- */
   /* make a movie */
   if (!strcmp(progname, "stalkmov")) {
      
#define NFRAMES 60
      for (i=0; i<NFRAMES; i++) {
         printf (" doing frame %d of %d \n", i, NFRAMES);
         tmp = ((double) i) / ((double) NFRAMES);
         re_center = -0.6;
         im_center = 0.0;
         width = 2.8;


         mandelbrot_stalk (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.3-1.5*tmp, 0.6*tmp); 
   

         /* dump the floating point data */
         sprintf (buff, "%s%d", argv[1], i);
         if ( (fp = Fopen (buff, ".flo")) == NULL) {
            printf (" File open failure for %s \n", buff);
            return 1;
         }
         fprintf (fp, "%d %d\n", data_width, data_height);
         fwrite (data, sizeof(float), data_width*data_height, fp);
         fclose (fp);
      }
      return 0;
   }
   
   /* ---------------------------------------------------- */
   /* make a movie */
   if (!strcmp(progname, "circmov")) {
      
      for (i=0; i<NFRAMES; i++) {
         printf (" doing frame %d of %d \n", i, NFRAMES);
         tmp = (((double) i) +0.25) / ((double) NFRAMES);
/*
         re_center = 0.0;
         im_center = 0.0;
         width = 1.0;
         width = 0.4; 
         circle_out (data, data_width, data_height,
                     re_center, im_center, width, itermax, tmp);
*/

         re_center = 0.5;
         im_center = 0.0;
         width = 1.0;
         circle_in (data, data_width, data_height,
                     re_center, im_center, width, itermax, tmp);

         /* dump the floating point data */
         sprintf (buff, "%s%d", argv[1], i);
         if ( (fp = Fopen (buff, ".flo")) == NULL) {
            printf (" File open failure for %s \n", buff);
            return 1;
         }
         fprintf (fp, "%d %d\n", data_width, data_height);
         fwrite (data, sizeof(float), data_width*data_height, fp);
         fclose (fp);
      }
      return 0;
   }
   
   /* dump the floating point data */
   if ( (fp = Fopen (argv[1], ".flo")) == NULL) {
      printf (" File open failure for %s.flo\n", argv[1]);
      return 1;
   }
   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (argv[1], ".txt")) == NULL) {
      printf (" File open failure for %s.txt\n", argv[1]);
      return 1;
   }
   fprintf (fp, "# Mandelbrot interior\n");
   fprintf (fp, "# floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "# width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "width %d height: %d \n", data_width, data_height);
   fprintf (fp, "height %d \n", data_height);
   fprintf (fp, "# Center, Real Part = %f \n", re_center);
   fprintf (fp, "# Center, Imaginary Part = %f \n", im_center);
   fprintf (fp, "# Width = %f \n", width);
   fprintf (fp, "# Data stored in natural units divided by %f\n", renorm);
   fprintf (fp, "# %d random samples taken per pixel\n", itermax);
   fprintf (fp, "# each sample iterated %d times\n", LOOP_COUNT);
   fprintf (fp, "# Therefore, each pixel normalized by %d * %d * %f\n",
		   itermax, LOOP_COUNT, renorm);
   fclose (fp);

   free (data);

   return 0;
}

/* --------------------------- END OF LIFE ------------------------- */
