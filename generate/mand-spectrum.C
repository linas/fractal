/*
 * mand-spectrum.C
 *
 * FUNCTION:
 * Explore spectral asymetry type sums of mandelbrot interior
 * Also includes interior/exterio decision routines.
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
/* This routine fills in the interior and exterior of the mandelbrot set 
 * using derivitivee w.r.t c (the infintessimal flow) to obtain values.
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
   double	re_c, im_c;
   double	re, im, tmp;
   double	dre, dim;
   double	ddre, ddim;
   double	d3re, d3im;
   double	d4re, d4im;
   double	dare, daim;
   double	dfre, dfim;
   double	dofre, dofim;
   int		loop;
   double modulus, phi, phip, phi3, frac;
   double escape_radius = 1000.0;
   double pot, ren, tl, otl;

   ren = log( log (escape_radius)) / log(2.0);
   tl = log(2.0);
   otl = 1.0/ log(2.0);
   
   /* adjust the iteration count to show the correct values */
   // itermax +=1;

   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = -1.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re_c = re_position;
         im_c = im_position;
         re = re_c;
         im = im_c;
         dre = 1.0;
         dim = 0.0;
         ddre = 0.0;
         ddim = 0.0;
         d3re = 0.0;
         d3im = 0.0;
         d4re = 0.0;
         d4im = 0.0;
         for (loop=1; loop <itermax; loop++) {
            /* compute fourth derivative */
            tmp = re*d4re - im*d4im + 4.0*(dre*d3re - dim*d3im);
            tmp += 3.0 * (ddre*ddre - ddim*ddim);
            tmp *= 2.0;
            d4im = re*d4im + im*d4re + 4.0*(dre*d3im + dim*d3re); 
            d4im += 6.0 * ddre * ddim;
            d4im *= 2.0;
            d4re = tmp;

            /* compute third derivative */
            tmp = 2.0 * (re*d3re - im*d3im + 3.0*(dre*ddre - dim*ddim));
            d3im = 2.0 * (re*d3im + im*d3re + 3.0*(dre*ddim + dim*ddre)); 
            d3re = tmp;

            /* compute second derivative */
            tmp = 2.0 * (re*ddre - im*ddim + dre*dre - dim*dim);
            ddim = 2.0 * (re*ddim + im*ddre + 2.0 * dre*dim);
            ddre = tmp;

            /* compute infinitessimal flow */
            tmp = 2.0 * (re*dre - im*dim) +1.0;
            dim = 2.0 * (re*dim + im*dre);
            dre = tmp;

            /* compute iterate */
            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         modulus = (re*re + im*im);
         modulus = sqrt (modulus);
         frac = log (log (modulus)) *otl;

         /* frac is the renormalized iteration count */
         frac = ((double) loop) - frac + 1.0; 

         /* pot is the duoady-hubbard potential */
         pot = exp (-tl*frac);

         /* compute d|z|/dc / |z| */
         dare = re*dre / (re*re+im*im);
         daim = im*dim / (re*re+im*im);

         /* compute dfrac/dc = d|z|/dc / |z| log|z| */
         dfre = tl*re*dre / ((re*re+im*im) *log (modulus));
         dfim = tl*im*dim / ((re*re+im*im) *log (modulus));

         /* compute 1/ (dfrac/dc)  */
         dofre = dfre / (dfre*dfre + dfim*dfim);
         dofim = -dfim / (dfre*dfre + dfim*dfim);

         /* compute zprime/z */
         tmp = re*dre + im*dim;   /* divergence */
         dim = re*dim - im*dre;   /* curl */
         dre = tmp;
         dre /= (re*re + im*im);
         dim /= (re*re + im*im);

         /* compute zprimeprime/z */
         tmp = re*ddre + im*ddim;   /* divergence */
         ddim = re*ddim - im*ddre;   /* curl */
         ddre = tmp;
         ddre /= (re*re + im*im);
         ddim /= (re*re + im*im);

         /* compute z'''/z */
         tmp = re*d3re + im*d3im;   /* divergence */
         d3im = re*d3im - im*d3re;   /* curl */
         d3re = tmp;
         d3re /= (re*re + im*im);
         d3im /= (re*re + im*im);

         /* compute z'(4)/z */
         tmp = re*d4re + im*d4im;   /* divergence */
         d4im = re*d4im - im*d4re;   /* curl */
         d4re = tmp;
         d4re /= (re*re + im*im);
         d4im /= (re*re + im*im);


         /* phase */
         /* just remember that gradient is contravarient
          * i.e. grad = 2 d-bar so flip sign of y 
          */
         phi = atan2 (-dim, dre);
         if (0.0 > phi) phi += 2.0*M_PI;
         phi /= 2.0*M_PI;

         phip = atan2 (-ddim, ddre);
         if (0.0 > phip) phip += 2.0*M_PI;
         phip /= 2.0*M_PI;

         phi3 = atan2 (-d3im, d3re);
         if (0.0 > phi3) phi3 += 2.0*M_PI;
         phi3 /= 2.0*M_PI;

         modulus = sqrt (dre*dre+dim*dim);
         modulus = sqrt (ddre*ddre+ddim*ddim);
         // modulus = sqrt (d4re*d4re+d4im*d4im);
         modulus /= (double) loop;
         modulus *= log((double) loop);
         modulus *= log((double) loop);

         modulus = sqrt (dofre*dofre+dofim*dofim);
         modulus = 1.0/sqrt (dre*dre+dim*dim);
         modulus = 1.0/sqrt (dare*dare+daim*daim);
        
         modulus = 0.25*(dre*dre+dim*dim)/(re*re+im*im);
         modulus /= log (sqrt(re*re+im*im));
         modulus /= log (sqrt(re*re+im*im));

         /* modulus of gradient of duoady-hubbard potential */
         modulus = pot;
         modulus /= log (sqrt(re*re+im*im));
         modulus *= sqrt(dre*dre+dim*dim);

if (loop>=itermax) phi = 0.0;
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine attempts to determine if 
 * a given point is inside or outside the mandelbrot set.
 */

void mandelbrot_decide (
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
   double	dre, dim;
   double	ddre, ddim;
   double	zppre, zppim;
   int		loop;
   double 	modulus, limit;
   double escape_radius = 50.0;
   int		state;
   double	*limits;
   clock_t	start, stop;
   int          hunds;

#define OUTSIDE    1
#define UNDECIDED  6
#define INSIDE     4
#define ERROR      8

   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = (double) UNDECIDED;

   limits = (double *)malloc ((itermax+1) * sizeof (double));
   for (i=1; i<itermax; i++) {
      double li = log ((double)i );
      limits[i] = 1.0;                            // too conservative
      limits[i] = (double)i;                      // too liberal
      limits[i] = sqrt ((double)i);               // too conservative
      limits[i] = ((double)i) / log ((double)i ); // too liberal
      limits[i] = ((double)i) / (li*li);          // just right !!
   }

   start = clock();
   hunds = 0;
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re = re_position;
         im = im_position;
         dre = 1.0;
         dim = 0.0;
         ddre = 0.0;
         ddim = 0.0;
         state = UNDECIDED;
         for (loop=1; loop <itermax; loop++) {

            /* compute second derivative */
            tmp = 2.0 * (re*ddre - im*ddim + dre*dre - dim*dim);
            ddim = 2.0 * (re*ddim + im*ddre + 2.0 * dre*dim);
            ddre = tmp;

            /* compute infinitessimal flow */
            tmp = 2.0 * (re*dre - im*dim) +1.0;
            dim = 2.0 * (re*dim + im*dre);
            dre = tmp;

            /* compute iterate */
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            modulus = (re*re + im*im);

            /* if point is outside the escape radius, then we've
             * determined that this is an exterior point.  If it
             * was ever marked as an interior point, that would be an
             * error.
             */
            if (modulus > escape_radius*escape_radius) {
               if ((UNDECIDED == state) || (OUTSIDE == state)) {
                  state = OUTSIDE;
               } else {
                  state = ERROR;
                  if (ERROR < loop) state = loop;
               }
               break; 
            }

            /* compute zprimeprime/z */
            zppre = re*ddre + im*ddim;   /* divergence */
            zppim = re*ddim - im*ddre;   /* curl */
            zppre /= (re*re + im*im);
            zppim /= (re*re + im*im);
   
            modulus = sqrt (zppre*zppre+zppim*zppim);
            // modulus /= (double) loop;
            limit = limits[loop];

            /* check to see if we can identify this as an interior point */
            if ((10 < loop) && (limit > modulus) && (UNDECIDED == state)) {
               state = INSIDE;
// record time saved ...
// state = loop;
break;
            }

         }    
        
         if (ERROR > state) {
            glob [i*sizex +j] = ((double) state) / 8.1;
         } else {
            glob [i*sizex +j] = ((double) state) / ((double)itermax);
         }

         re_position += delta;

         /* Measure cpu time spent.  On Linux, there are a million
          * clocks per sec, which will overflow in about an hour.
          * so we have to measure time frequently!
          */
         if (j%100==0) {
            stop = clock();
            hunds += (stop-start) / (CLOCKS_PER_SEC/100);
            start = stop;
         }
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
   stop = clock();
   hunds += (stop-start) / (CLOCKS_PER_SEC/100);
   start = stop;

   free (limits);

   /* do some counting. */
   { 
   int inside=0, uncertain=0;
   double rem;
   for (i=0; i<globlen; i++) {
      if ((((double) UNDECIDED)+0.1 > 8.1*glob[i]) &&
          (((double) UNDECIDED)-0.1 < 8.1*glob[i])) { uncertain ++; }
      if ((((double) INSIDE)+0.1 > 8.1*glob[i]) &&
          (((double) INSIDE)-0.1 < 8.1*glob[i])) {inside ++; }
   }
   rem = ((double) uncertain) / ((double)(uncertain+inside));
   printf (">> %d %d %d %d %f %d\n", 
       itermax, uncertain, inside, 
       uncertain+inside, rem, hunds);
   }
}

/* --------------------------- END OF LIFE ------------------------- */
