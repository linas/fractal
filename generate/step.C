/*
 * step.C
 *
 * FUNCTION:
 * debug
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 */

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


main () 
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
   double phi=0.0, phi_last, phi_c, h_phi_c;
   int wind =0;

   int itermax =200;

         re_c = -0.2;
         // im_c = 0.842;
         im_c = 0.843;

         phi_c = atan2 (im_c, re_c);
         if (0.0 > phi_c) phi_c += 2.0*M_PI;
         h_phi_c = 0.5*phi_c;

         re = 0.0;
         im = 0.0;
         re_prev = 0.0;
			im_prev = 0.0;
         phi_last = -0.01;
         wind = 0;
printf ("\n phi_c= %12.8g   c=(%g %g)\n", phi_c, re_c, im_c);
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

				// The below cancels out horns 
            // if (phi < phi_last) wind +=2;

printf ("n=%d %12.8g %12.8g %12.8g %12.8g w=%d\n", 
loop, phi, phi_last, 2.0*phi_last-phi, 
phi+((double)wind)*M_PI, wind);

				/* The following is critical for going from discontinuous
				 * power-two bars to smooth func. */
            /* if northern half else southern half */
            if (M_PI > phi_c) { 
               while (M_PI < phi) {
                  phi -= 2*M_PI;
                  wind +=2;
               }
            } else {
               while (M_PI > phi) {
                  phi += M_PI;
                  wind --;
               }
            }

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

}

/* --------------------------- END OF LIFE ------------------------- */
