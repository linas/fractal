/*
 * smooth.c
 *
 * FUNCTION:
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 */

#include <stdio.h>
#include <math.h>

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_out (void)
{
   int	i;
   double corr;
   
   double		re_position, im_position;
   double		re, im, tmp;
   int		loop;
   double modulus, ln, lnp1, lnp2;
   double escape_radius = 33.5;
   double fix, tl;
   double prev = 0.0;
   int nprev = 0;
   double alpha, beta;
   double gamma, delta;
   double da;
   double rbar, rbar1;
   double mu;

   tl = 1.0/ log(2.0);
   
   
   re_position = -1.0;
   im_position = 0.5;

   for (i=0; i<500; i++) {
      escape_radius = 2.0 + 0.05 * (double) i*i*i;
      fix = log( log (escape_radius)) / log(2.0);

      re = re_position;
      im = im_position;
      for (loop=1; loop <800; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         modulus = (re*re + im*im);
         if (modulus > escape_radius*escape_radius) break;
      }    

      modulus = sqrt (modulus);
      ln = log (log (modulus)) *tl;

      mu = (double)loop - ln;

      tmp = re*re - im*im + re_position;
      im = 2.0*re*im + im_position;
      re = tmp;
      modulus = (re*re + im*im);
      modulus = sqrt (modulus);
      lnp1 = log (log (modulus)) *tl;

      rbar = lnp1 -ln -1.0;

      tmp = re*re - im*im + re_position;
      im = 2.0*re*im + im_position;
      re = tmp;
      modulus = (re*re + im*im);
      modulus = sqrt (modulus);
      lnp2 = log (log (modulus)) *tl;

      rbar1 = lnp2 -lnp1 -1.0;

      /* second order correction */
      alpha = ln-fix;
      beta = 1.0 - alpha;
      gamma = 4.0 * alpha*beta;
      delta = 1.0 - gamma;

      da = delta * alpha;
      corr = gamma *ln + delta*fix;
      corr = mu - rbar - rbar1;
       
      nprev = loop - nprev;
      prev = corr -prev; 
      printf ("r=%f n=%d mu=%f Ln=%f r=%e r1=%e corr=%f diff=%e\n",
          escape_radius, loop, mu, ln, rbar, rbar1, corr, prev);
      nprev = loop;
      prev = corr;
   }
}

void main (int argc, char *argv[]) 
{
   
/*
   if (5 > argc) {
      fprintf (stderr, "Usage: %s \n", argv[0]);
      exit (1);
   }
*/


   mandelbrot_out ();
}

/* --------------------------- END OF LIFE ------------------------- */
