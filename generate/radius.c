/*
 * radius.c
 *
 * FUNCTION:
 * Find bud radius.  Itdoes this by doing iterations along a set of radial
 * lines coming out from a center.  When iterated sufficiently, the iterator
 * either settles down to a cycle (in which case the pont is inside) or it
 * escapes (i.e. we have an over-estimate for the radius).
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 */

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*-------------------------------------------------------------------*/
/* this routine measures a bud radius using */
/* the classic algorithm */

void measure_radius (
   double  	*glob,
   unsigned int sizea,
   unsigned int sizer,
   double	re_center,
   double	im_center,
   double	rmin,
   double	rmax,
   double	epsilon,
   int		itermax)
{
   unsigned int	i,j;
   double	deltar, rad, theta=0.0;
   double	radius_cand;   /* minimum possible radius */
   double	radius_outer;  /* escape at this radius */
   double	re_cg, im_cg;
   double	re_outer_cg, im_outer_cg;
   double	radius_avg=0.0;     /* average over all samples */
   double	radius_min, radius_max;  /* min and max radius over entire bud */
   double	radius_outer_avg=0.0;
   double	radius_outer_min, radius_outer_max;
   double	si, co;
   double	re_last, im_last;
   double	re_c, im_c;
   double	re, im, tmp;
   int		loop=0, loop_cand=0;
   double 	modulus=0.0;
   double	escape_radius = 3.1;
   double	esq;
   double	dist=0.0;
   double	closest=0.0;

   esq = epsilon*epsilon;

   deltar = (rmax-rmin) / (double) sizer;
   
   for (i=0; i<sizea; i++) glob [i] = 0.0;
   
   re_cg = im_cg = 0.0;
   radius_avg = 0.0;
   radius_min = 1e30;
   radius_max = -1e30;

   re_outer_cg = im_outer_cg = 0.0;
   radius_outer = 0.0;
   radius_outer_avg = 0.0;
   radius_outer_min = 1e30;
   radius_outer_max = -1e30;

   for (i=0; i<sizea; i++) {
      theta = ((double) i) / ((double) sizea);
      theta *= 2.0 * M_PI;
      si = sin (theta);
      co = cos (theta);

      rad = rmin;
      radius_cand = rad;
      for (j=0; j<sizer; j++) {
         re_c = re_center + rad * co;
         im_c = im_center + rad * si;
   
         closest = 1.0e30;
         re = re_c;
         im = im_c;
         for (loop=1; loop <itermax; loop++) {
            re_last = re;
            im_last = im;

            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;

            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;

            tmp = re*re - im*im + re_c;
            im = 2.0*re*im + im_c;
            re = tmp;

            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) {
               radius_outer = rad; /* radius must be smaller than this */
               j+= sizer; /* break middle loop */
               break;
            }

            dist = (re-re_last)*(re-re_last) + (im-im_last)*(im-im_last);
            if (dist < closest) closest = dist;
            if (dist < esq) {
               radius_cand = rad;
               loop_cand = loop;
               break;
            }
         }    
         printf ("%14.10g	%14.10g	%d\n", rad, closest, loop);
   
         rad += deltar;
      }
      re_cg += radius_cand * co;
      im_cg += radius_cand * si;
      radius_avg += radius_cand;
      if (radius_min > radius_cand) { radius_min = radius_cand; }
      if (radius_max < radius_cand) { radius_max = radius_cand; }
      glob [i] = radius_cand; 

      re_outer_cg += radius_outer * co;
      im_outer_cg += radius_outer * si;
      radius_outer_avg += radius_outer;
      if (radius_outer_min > radius_outer) { radius_outer_min = radius_outer; }
      if (radius_outer_max < radius_outer) { radius_outer_max = radius_outer; }

      printf ("%d	%14.10f	%14.10f	%14.10f %g	%d\n", 
            i, theta, radius_cand, radius_outer, radius_outer-radius_cand, loop_cand); 
   }
   radius_avg /= (double) sizea;
   re_cg /= (double) sizea;
   im_cg /= (double) sizea;
   
   radius_outer_avg /= (double) sizea;
   re_outer_cg /= (double) sizea;
   im_outer_cg /= (double) sizea;
   
   printf ("# ravg = %14.10f center o gravity = ( %14.10f %14.10f )\n", 
        radius_avg, re_cg, im_cg);
   printf ("# rout = %14.10f outer center o g = ( %14.10f %14.10f )\n", 
        radius_outer_avg, re_outer_cg, im_outer_cg);
   printf ("# center= %14.10f %14.10f diam = %14.10f \n", 
        re_center+re_cg, im_center+im_cg, 2.0*radius_avg);
   printf ("# outcen= %14.10f %14.10f out diam = %14.10f \n", 
        re_center+re_outer_cg, im_center+im_outer_cg, 2.0*radius_outer_avg);
   printf ("# radius min, max= %14.10f %14.10f half-diff=%14.10f \n", 
        radius_min, radius_max, 0.5*(radius_max-radius_min));
   printf ("# outer  min, max= %14.10f %14.10f outer-half=%14.10f \n", 
        radius_outer_min, radius_outer_max, 0.5*(radius_outer_max-radius_outer_min));
}


/*-------------------------------------------------------------------*/

int 
main (int argc, char *argv[]) 
{
   double	*data;		/* my data array */
   unsigned int	nphi, nr;
   double	re_center, im_center, rmin, rmax;
   int		itermax;
   double	epsilon;
   
   if (5 > argc) {
      fprintf (stderr, "Usage: %s <n_phi> <n_r> <niter> <epsilon> [<centerx> <centery> <rmin> <rmax>]\n", argv[0]);
      exit (1);
   }

   itermax = 1;

   nphi = atoi (argv[1]);
   nr = atoi (argv[2]);
   itermax = atoi (argv[3]);
   epsilon = atof (argv[4]);

   data = (double *) malloc (nphi*sizeof (double));

   /* do bud at top, the 3-loop */
   re_center = -0.125;
   im_center = 0.7445;
   rmin = 0.09;
   rmax = 0.1;
   if (argc >= 8) {
      re_center = atof (argv[5]);
      im_center = atof (argv[6]);
      rmin = atof (argv[7]);
      rmax = atof (argv[8]);
   }

   printf ("# \n");
   printf ("# measurement of loop=3 radius \n");
   printf ("# \n");
   printf ("# nphi=%d nr=%d iter=%d eps=%g cent=(%14.10f %14.10f) rmin=%f rmax=%f\n", 
        nphi, nr, itermax, epsilon, re_center, im_center, rmin, rmax);
   printf ("# \n");
   printf ("#i	theta		radius		radius outer	outer-inner	loop count\n"); 
   printf ("# \n");

   measure_radius (data, nphi, nr, re_center, im_center,
            rmin, rmax, epsilon, itermax);
   

   free (data);

   return 0;
}

/* --------------------------- END OF LIFE ------------------------- */
