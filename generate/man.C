/* 
 * NAME:
 * man.c
 *
 * FUNCTION:
 * generate pixmaps showing Poincare recurrence time for mandelbrot set
 *
 * HISTORY:
 * created Linas Vepstas Jan 1991
 * spruced up Linas Vepstas July 1993
 * fixed major performance bug (AIX fmod() stinks!) -- January 1994
 * added logistic map -- February 1994
 * added classic mandelbrot -- June 1995
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "image.h"
#include "util.h"

/*-------------------------------------------------------------------*/

#define EPSILON  	0.0003
/* 
#define SETTLE_TIME 	500
*/
#define SETTLE_TIME 	800

/*-------------------------------------------------------------------*/

typedef struct _circ_data {
   double epsilon;
   double normalization;
   int settle_time;
   int itermax;
} CircleData;

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the classic z**2+c
 * mandelbrot map
 */

#define MNORM 5.0
double mandel_cross = 1.0;

float mandelbrot_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y, tmp;
   double	cx, cy;
   double	dist, dx, dy;
   double	xpoint, ypoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   cx = K * cos (omega) - 1.0;
   cy = K * sin (omega);

/*
   theta = omega / M_PI - 1.0;
   r = K / (K*K + theta*theta);
   phi = theta / (K*K+theta*theta);

   cx = r * cos (M_PI*phi);
   cy = r * sin (M_PI*phi);
*/

   x = cx;
   y = cy;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      tmp = x*x - y*y + cx;
      y = mandel_cross * 2.0 * x * y + cy;
      x = tmp;
      if (10.0 < (x*x+y*y)) return 0.0;
   }


   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   ypoint = y;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      tmp = x*x - y*y + cx;
      y = mandel_cross * 2.0 * x * y + cy;
      x = tmp;

      dx = x-xpoint;
      dy = y-ypoint;
      dist = dx*dx + dy*dy;
      if (dist < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (MNORM * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
   
void do_mandel (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;
   char buff[100];
   int iii;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   /*   The BIG picture
   */
   omega_min = 0.0;
   omega_max = 2.0 * M_PI;
   K_min = 2.0; 
   K_max = 0.0;

/*
   omega_min = 0.8 * M_PI;
   omega_max = 1.2 * M_PI;
   K_min = 1.9; 
   K_max = 1.3;
*/
   
   dat.itermax = 30;
/*
mandel_cross = +0.1;
for (iii=0; iii<6; iii++) {
mandel_cross -= 0.1;
*/
{
mandel_cross = 1.0; iii=3;

   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              (DensityCB) mandelbrot_poincare_recurrance_time, &dat);

   
   /* dump the floating point data */
   sprintf (buff, "%s%d", filename, iii);
   if ( (fp = Fopen (buff, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (buff, ".txt")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "mangled mandelbrot map\n");
   fprintf (fp, "x(n+1) = x(n)*x(n) - y(n)*y(n) + ReC \n");
   fprintf (fp, "y(n+1) = 2*i*x(n)*y(n) * D + ImC \n");
   fprintf (fp, "where D = %f \n", mandel_cross);
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", MNORM);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "lambda on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);
}
}

/* ----------------------- end of file  ---------------------- */
   
