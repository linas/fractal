
/*
 * plot ray vs. theta
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
ray (double *bin, int nbin)
{
   double re, im;
   double dre, dim;
   double phi;
   double escape_radius, esq;
   int i, j, loop, itermax;
   int nx, ny;
   double re_c, im_c;
   double delta_re, delta_im;
   double tmp, modulus;
   int *cnt;
   int ibin;
   double u,v, r;

   nx = ny = 800;
   itermax = 100000;
   escape_radius = 1e10;

   delta_re = 3.0/((double) nx);
   delta_im = 3.0/((double) ny);

   esq = escape_radius*escape_radius;

   // malloc a local counter array
   cnt = (int *) malloc (nbin * sizeof (int));

   // initialize the bins
   for (i=0; i<nbin; i++){
      bin[i] = 0.0;
      cnt[i] = -1;
   }

   // walk over pixel grid
   re_c = -2.1;
   for (i=0; i<nx; i++) 
      {

      if (0 == i%10) printf ("# doing row %d\n", i);

      im_c = -1.5;
      for (j=0; j<ny; j++) 
         {

         // if we are inside the cardiod, then don't even bother
         // this saves a significant amount of cpu
         re = 0.25 - re_c;
         im = -im_c;
         u = sqrt (0.5*(re + sqrt (re*re + im*im)));
         v = 0.5 * im / u;
         u = 0.5 - u;
         r = sqrt (u*u + v*v);
         if (0.5 > r) {
            im_c += delta_im;
            continue;
         }


         // initialize the iterators
         re = 0.0;
         im = 0.0;
         dre = 0.0;
         dim = 0.0;

         // iterate until escape
         for (loop=0; loop<itermax; loop++) {
            tmp = 2.0 * (re*dre - im*dim) + 1.0;
            dim = 2.0 * (re*dim + im*dre);
            dre = tmp;

            tmp = re*re - im*im + re_c;
            im = 2.0 * re*im + im_c;
            re = tmp;

            modulus = re*re + im*im;
            if (modulus > esq) break;
         }

         // do the math only is this point escaped
         if (modulus > esq) {
            phi = atan2 (im, re);
            phi /= 2.0 * M_PI;
            phi += 0.5;
            phi *= nbin;
            ibin = (int) phi;
            if (ibin == nbin) ibin = 0;

            // update only if we got closer to the edge
            if (loop > cnt[ibin]) {
               cnt[ibin] = loop;
               bin[ibin] = atan2 (dim, dre) / (2.0*M_PI) + 0.5;
            }
         }

         im_c += delta_im;
      }
      re_c += delta_re;
   }

   free (cnt);
}


int
main (int argc, char *argv[])
{

   double *map;
   int i, nbin;
   double phi;

   nbin = 400;

   map = (double *) malloc (nbin * sizeof(double));

   ray (map, nbin);

   for (i=0; i<nbin; i++) {
      phi = ((double) i) / ((double) nbin);
      printf ("%d	%f	%f\n", i, phi, map[i]);
   }

   return 0;
}
