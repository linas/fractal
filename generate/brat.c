/*
 * brat.c
 *
 * FUNCTION:
 * Explore Hausdorf measure of mandelbrot set.
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
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_out (
   float  	*glob,
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   unsigned int	i,j, globlen;
   double		re_start, im_start, delta;
   double		re_position, im_position;
   double		re, im, tmp;
   int		loop;
   double modulus=0.0, frac;
   double escape_radius = 3.1;
   double ren, tl;

   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0/ log(2.0);
   

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
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         modulus = sqrt (modulus);
         frac = log (log (modulus)) *tl;

         /* frac =  Re (c/z*z) */
         tmp = (re*re-im*im);
         im = 2.0*re*im;
         re = tmp;
         frac = re*re_position + im*im_position;
         frac /= re*re+im*im;
         if (0.0 > frac) frac = -frac;
         
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
         glob [i*sizex +j] = ((double) loop); 
/*
         glob [i*sizex +j] = ((float) (loop%10)) / 10.0; 
if (loop == itermax) {
glob[i*sizex+j] = 0.0; } else {glob[i*sizex+j]=0.9999;}
*/

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the mandelbrot set using 
 * the classic algorithm. The derivative (infinitessimal flow) 
 * is used to play games.
 */

void dmandelbrot_out (
   float  	*glob,
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   unsigned int	i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re, im, tmp;
   double	dre, dim;
   double	ddre, ddim;
   double	d3re, d3im;
   double	d4re, d4im;
   int		loop;
   double modulus, phi, phip, phi3, frac;
   double escape_radius = 3450.0;
   double ren, tl;

   ren = log( log (escape_radius)) / log(2.0);
   tl = 1.0/ log(2.0);
   

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
         re = re_position;
         im = im_position;
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
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            modulus = (re*re + im*im);
            if (modulus > escape_radius*escape_radius) break;
         }    

         modulus = (re*re + im*im);
         modulus = sqrt (modulus);
         frac = log (log (modulus)) *tl;

         frac = ((double) loop) - frac + 1.0; 

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

         /* compute z(4)/z */
         tmp = re*d4re + im*d4im;   /* divergence */
         d4im = re*d4im - im*d4re;   /* curl */
         d4re = tmp;
         d4re /= (re*re + im*im);
         d4im /= (re*re + im*im);


         /* phase */
         phi = atan2 (dim, dre);
         phi += M_PI;
         phi /= 2.0*M_PI;

         phip = atan2 (ddim, ddre);
         phip += M_PI;
         phip /= 2.0*M_PI;

         phi3 = atan2 (d3im, d3re);
         phi3 += M_PI;
         phi3 /= 2.0*M_PI;

         modulus = sqrt (dre*dre+dim*dim);
         modulus = sqrt (ddre*ddre+ddim*ddim);
         // modulus = sqrt (d4re*d4re+d4im*d4im);
         modulus /= (double) loop;
         modulus *= log((double) loop);
         modulus *= log((double) loop);

        
         glob [i*sizex +j] = modulus;

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
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   unsigned int	i,j, globlen;
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
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
      stop = clock();
      hunds += (stop-start) / 10000;
      start = stop;
   }
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

/*-------------------------------------------------------------------*/
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_inverse (
   float  	*glob,
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   unsigned int	i,j, globlen;
   double		re_start, im_start, delta;
   double		re_position, im_position;
   double		re, im, tmp;
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

         ire = sizex * (re + 0.5);
         if (0 > ire) ire = 0;
         if (sizex <= ire) ire = sizex-1;

         im *= ((double) sizey)/ ((double) sizex);
         jim = sizey * (im + 0.5);

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
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   unsigned int	i,j, globlen;
   double		re_start, im_start, delta;
   double		re_position, im_position;
   double		re, im, tmp;
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
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double	stalkx,
   double	stalky)
{
   unsigned int	i,j, globlen;
   double		re_start, im_start, delta;
   double		re_position, im_position;
   double		re, im, tmp, tmpx, tmpy;
   double		visited_x, visited_y;
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
#endif MAXDIST
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
#endif CIRCSTALK

           
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
#endif MAXDIST

           


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
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	time)
{
   unsigned int	i,j, globlen;
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
   unsigned int sizex,
   unsigned int sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax,
   double 	time)
{
   unsigned int	i, globlen;
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
   isamp = CBOX_IM_SLOPE*CBOX_RE_SLOPE / (x_width * y_width);
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

#define BBOX_IM_SLOPE 3.0
#define BBOX_IM_CEPT -1.5
#define BBOX_RE_SLOPE 3.0
#define BBOX_RE_CEPT -2.1

void mandelbrot_measure (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_offset (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */
{
   unsigned int	globlen;
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
   
   isamp = ((double) BBOX_IM_SLOPE*BBOX_RE_SLOPE) / (x_width * y_width);
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

void mandelbrot_age (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_lyapunov (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_migrate (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_squige (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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
            last = ++ count [vert_pix*sizex + horiz_pix];
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

void mandelbrot_phase (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_orig (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void mandelbrot_next (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void whack (glob, sizex, sizey,
               re_center, im_center, width, itermax, renorm)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
double		renorm;
{
   unsigned int	globlen;
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
   
   isamp = BBOX_IM_SLOPE*BBOX_RE_SLOPE / (x_width * y_width);
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

void random_out (glob, sizex, sizey,
               re_center, im_center, width, itermax)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
/* This routine fills in the exterior of the mandelbrot set using 
 * an algorithm which is, well, a bit differnet. 
 * I suspect a problem with the original algorithm is that the seed
 * pixels always start on grid boundaries.  The algorithm below
 * provides a random seed instead, and assigns that to the pixel.
 * Like the traditional exterior coloring algorithms, the number of
 * iterations a pixel goes through before it escapes is counted.
 */
{
   unsigned int	globlen;
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

extern FILE *Fopen();

int 
main (int argc, char *argv[]) 
{
   float	*data;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   double	re_center, im_center, width;
   int		itermax;
   double	renorm, tmp;
   FILE		*fp;
   int		i;
   char buff [80];
   
   if (5 > argc) {
      fprintf (stderr, "Usage: %s <filename> <width> <height> <niter> [<centerx> <centery> <size>]\n", argv[0]);
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

   if (8 == argc) {
      re_center = atof (argv[5]);
      im_center = atof (argv[6]);
      width = atof (argv[7]);
   }

   /* Do the interior now */
   renorm = 1.0;

   if (!strcmp(argv[0], "brat"))
   mandelbrot_out (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(argv[0], "inf"))
   dmandelbrot_out (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(argv[0], "decide"))
   mandelbrot_decide (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(argv[0], "manvert"))
   mandelbrot_inverse (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(argv[0], "measure"))
   mandelbrot_measure (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(argv[0], "offset"))
   mandelbrot_offset (data, data_width, data_height,
                  0.5, 0.0, 3.0, itermax, renorm); 
   
   if (!strcmp(argv[0], "mstop"))
   mandelbrot_stop (data, data_width, data_height,
                  re_center, im_center, width, itermax); 
   
   if (!strcmp(argv[0], "stalk"))
   mandelbrot_stalk (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.0, 0.0); 
   
   if (!strcmp(argv[0], "orig"))
   mandelbrot_orig (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(argv[0], "phase"))
   mandelbrot_phase (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(argv[0], "age"))
   mandelbrot_age (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
#ifdef STRAIGHT
   if (!strcmp(argv[0], "lyapunov"))
   mandelbrot_lyapunov (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
#else
   if (!strcmp(argv[0], "lyapunov"))
   mandelbrot_lyapunov (data, data_width, data_height,
                  0.5, 0.0, 3.0, itermax, renorm); 
#endif
   

   if (!strcmp(argv[0], "migrate"))
   mandelbrot_migrate (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(argv[0], "squige"))
   mandelbrot_squige (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 

   if (!strcmp(argv[0], "next"))
   mandelbrot_next (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(argv[0], "whack"))
   whack (data, data_width, data_height,
                  re_center, im_center, width, itermax, renorm); 
   
   if (!strcmp(argv[0], "circout")) {
      re_center = 0.0;
      im_center = 0.0;
      width = 1.0;
      
      circle_out (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.33); 
   }
   
   if (!strcmp(argv[0], "circin")) {
      re_center = 0.5;
      im_center = 0.0;
      width = 1.0;
      
      circle_in (data, data_width, data_height,
                  re_center, im_center, width, itermax, 0.05); 
   }
   
   /* ---------------------------------------------------- */
   /* make a movie */
   if (!strcmp(argv[0], "stalkmov")) {
      
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
   if (!strcmp(argv[0], "circmov")) {
      
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
