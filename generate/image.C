/* 
 * NAME:
 * image.c
 *
 * FUNCTION:
 * generate pixmaps showing Poincare recurrence time for circle map.
 *
 * HISTORY:
 * created Linas Vepstas Oct 1989
 * major revs Linas Vepstas Jan 1991
 * spruced up Linas Vepstas July 1993
 * fixed major performance bug (AIX fmod() stinks!) -- January 1994
 * added logistic map -- February 1994
 * added classic mandelbrot -- June 1995
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "util.h"

extern void do_mandel (char * filename, int width, int height);
extern void do_cf_gap (char * filename, int width, int height);

typedef float (*DensityCB) (void*, double, double);

int param;

/*-------------------------------------------------------------------*/
/* Test routine, draws very pretty moire pattern */

void moire_parabola (float *glob, unsigned int sizex, unsigned int sizey)
{
   unsigned int	i,j, globlen;
   int		icen, jcen, ip, jp;
   float	globby;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   icen = (int) (0.65* sizex);
   jcen =  (int) (0.8 * sizey);
   for (i=0; i<sizey; i++)
      {
      for (j=0; j<sizex; j++) 
         {
         ip = i - icen; 
         jp = j - jcen;
         globby = (float) ((param*(ip*ip+jp*jp)) % 235);
         globby = globby / 235.0;
         glob [i*sizex +j] = globby; 
         }
      }
}

void do_moire (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   FILE		*fp;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   
   /* fill it in */
   moire_parabola (data, data_width, data_height);

   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "moire parabola\n");
   fclose (fp);
}

/*-------------------------------------------------------------------*/
/* 
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 */

void walk_rect (float *glob, 
                unsigned int sizex,  /* width, in pixels */
                unsigned int sizey,  /* height, in pixels */
                double x_min,        /* left side of pixmap */
                double x_max,        /* right side of pixmap */
                double y_min,        /* bottom of pixmap */
                double y_max,        /* top of pixmap */
                float (*callback)(void *, double, double),  /* callback */
                void * calldata)     /* static data */
{
   unsigned int	i,j, globlen;
   double	delta_x, delta_y;
   double x, y;
   
   
   /* initialize pixmap to zero */
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   /* compute slopes */
   delta_x = (x_max - x_min) / ((float) sizex);
   delta_y = (y_max - y_min) / ((float) sizey);
   
   /* loop over columns, then rows */
   for (i=0; i<sizey; i++) { 
      printf (" starting %d out of total %d \n", i, sizey);
      fflush (stdout);
      y = (((double) i) +0.5) * delta_y + y_min ;

      for (j=0; j<sizex; j++) {
         x = (((double) j) +0.5) * delta_x + x_min;

         /* get the values */
         glob [i*sizex +j] = (*(callback)) (calldata, x, y);
   
      }
   }
}
   
/*-------------------------------------------------------------------*/
/* 
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 * instead of filling in in a linear fashion, its done with a
 * trianglelar distortion.
 */

void walk_tri (float *glob, 
               unsigned int sizex,  /* width, in pixels */
               unsigned int sizey,  /* height, in pixels */
               double x_min,        /* left side of pixmap */
               double x_max,        /* right side of pixmap */
               double y_min,        /* bottom of pixmap */
               double y_max,        /* top of pixmap */
               float (*callback)(void *, double, double),  /* callback */
               void * calldata)     /* static data */
{
   unsigned int	i,j, globlen;
   double	delta_x, delta_y;
   double x, y;
   double x_start;
   
   
   /* initialize pixmap to zero */
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   /* compute slopes */
   delta_y = (y_max - y_min) / ((float) sizey);
   
   /* loop over columns, then rows */
   for (i=0; i<sizey; i++) { 
      printf (" starting %d out of total %d \n", i, sizey);
      y = (((double) i) +0.5) * delta_y + y_min ;

      delta_x = (x_max - x_min) / ((float) sizex);
      delta_x *= (y_max - y) / (y_max - y_min);

      x_start = x_min;
      x_start += ((y-y_min)/(y_max-y_min)) * (x_max-x_min) * 0.5;

      for (j=0; j<sizex; j++) {
         x = (((double) j) +0.5) * delta_x + x_start;

         /* get the values */
         glob [i*sizex +j] = (*(callback)) (calldata, x, y);
   
      }
   }
}
   
/*-------------------------------------------------------------------*/
/* 
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 * instead of filling in in a linear fashion, its done with a
 * trianglelar distortion.
 */

void walk_utri (float *glob, 
                unsigned int sizex,  /* width, in pixels */
                unsigned int sizey,  /* height, in pixels */
                double x_min,        /* left side of pixmap */
                double x_max,        /* right side of pixmap */
                double y_min,        /* bottom of pixmap */
                double y_max,        /* top of pixmap */
                float (*callback)(void *, double, double),  /* callback */
                void * calldata)     /* static data */
{
   unsigned int	i,j, globlen;
   double	delta_x, delta_y;
   double x, y;
   double x_start;
   
   
   /* initialize pixmap to zero */
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   /* compute slopes */
   delta_y = (y_max - y_min) / ((float) sizey);
   
   /* loop over columns, then rows */
   for (i=0; i<sizey; i++) { 
      printf (" starting %d out of total %d \n", i, sizey);
      y = (((double) i) +0.5) * delta_y + y_min ;

      delta_x = (x_max - x_min) / ((float) sizex);
      delta_x *= (y - y_min) / (y_max - y_min);

      x_start = x_min;
      x_start += ((y_max-y)/(y_max-y_min)) * (x_max-x_min) * 0.5;

      for (j=0; j<sizex; j++) {
         x = (((double) j) +0.5) * delta_x + x_start;

         /* get the values */
         glob [i*sizex +j] = (*(callback)) (calldata, x, y);
   
      }
   }
}
   
/*-------------------------------------------------------------------*/

typedef struct _circ_data {
   double epsilon;
   double normalization;
   int settle_time;
   int itermax;
} CircleData;


/*-------------------------------------------------------------------*/
/* 
 * This routine computes average winding number taken by
 * circle map iterator.
 */

float circle_winding_number (CircleData *dat,
                       double omega,
                       double K)

{
   double	x=0.0;
   int		iter;
   
  
   /* OK, now start iterating the circle map */
   for (iter=0; iter < dat->itermax; iter++) {
      x += omega - K * sin (2.0 * M_PI * x);
   }

   /* x is the normalized winding number */
   x /= (double) (dat->itermax);

/* ============
   if (x<0.0) x = 0.0;
   if (x>1.0) x = 1.0;
   if ((x<0.0) | (x>1.0)) {
      printf ("out of bounds omega=%f and K=%f W=%f \n", omega, K, x);
   }
   printf ("for omega=%f and K=%f W=%f \n", omega, K, x);
   ==============*/

   return ((float) x);
}
   
/*-------------------------------------------------------------------*/

#define EPSILON  	0.003
/* 
#define SETTLE_TIME 	500
*/
#define SETTLE_TIME 	800
#define NORMALIZATION 	400.0 

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 */

float circle_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	xpoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      x += omega - K * sin (2.0 * M_PI * x);
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      x += omega - K * sin (2.0 * M_PI * x);
      y = fabs (x-xpoint);
      y -= floor (y);
      if (y < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (NORMALIZATION * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 * for an imaginary coupling constant.
 */

#define INOR 10.0

float im_circle_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	tx, ty, ey;
   double 	cw, sw;
   double	xpoint, ypoint;
   double	dx, dy, d;
   int		iter;
   long		num_recurs, time_recur=0;

   cw = cos (2 * M_PI * omega);
   sw = sin (2 * M_PI * omega);
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 1.0;
   y = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      ey = exp (-K*y);
      tx = ey * (x * cw - y * sw);
      ty = ey * (x * sw + y * cw);
      x = tx;
      y = ty;
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   ypoint = y;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      ey = exp (-K*y);
      tx = ey * (x * cw - y * sw);
      ty = ey * (x * sw + y * cw);
      x = tx;
      y = ty;
      dx = x - xpoint;
      dy = y - ypoint;
      d = dx*dx + dy*dy;

      if (d < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (INOR * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 * for an complex coupling constant.
 */

double complex_phase =0.0;

#define CNOR 30.0
#define PHI complex_phase

float co_circle_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double ka, kb;
   double 	cw, sw;
   double	x, y;
   double	tx, ty, sy, cy, ey;
   double	xpoint, ypoint;
   double	dx, dy, d;
   int		iter;
   long		num_recurs, time_recur=0;

   cw = cos (2 * M_PI * omega);
   sw = sin (2 * M_PI * omega);

   ka = K * cos (PHI);
   kb = K * sin (PHI);
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 1.0;
   y = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      ey = exp (-kb*y);
      sy = sin (ka*y);
      cy = cos (ka*y);
      tx = (x * cw - y * sw);
      ty = (x * sw + y * cw);
      x = (tx * cy - ty * sy) * ey;
      y = (tx * sy + ty * cy) * ey;
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   ypoint = y;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      ey = exp (-kb*y);
      sy = sin (ka*y);
      cy = cos (ka*y);
      tx = (x * cw - y * sw);
      ty = (x * sw + y * cw);
      x = (tx * cy - ty * sy) * ey;
      y = (tx * sy + ty * cy) * ey;
      dx = x - xpoint;
      dy = y - ypoint;
      d = dx*dx + dy*dy;

      if (d < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (CNOR * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the logistic map
 * (feigenbaum map).
 */
#define NORMLOGIS 200.0

float logistic_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	xpoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = omega;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      x = K * x * (1.0 - x);
   }

   if (x > 100.0) return 0.0;
   if (x < -100.0) return 0.0;

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      x = K * x * (1.0 - x);

      y = fabs (x-xpoint);
      if (y < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (NORMLOGIS * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the logistic map
 * (feigenbaum map), clamped.
 */

float clamp_logistic_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	xpoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      x -= omega;
      /* x = K * x * (1.0 - x); */
      /* x = K * sin (M_PI*x);  */
      x = K * sin (M_PI*x) / (M_PI* x);   /* sinc function */
      /* x = K * exp (-x*x); */
      x -= floor (x);
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      x -= omega;
      /* x = K * x * (1.0 - x); */
      /* x = K * sin (M_PI*x); */
      x = K * sin (M_PI*x) / (M_PI* x);
      /* x = K * exp (-x*x); */
      x -= floor (x);

      y = fabs (x-xpoint);
      y -= floor (y);
      if (y < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (NORMLOGIS * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the logistic map
 * (feigenbaum map), clamped.
 */

float squash_logistic_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	xpoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* do this only  ... */
/*
   tmp = omega * K;
   omega = K - omega;
   K = tmp;
*/
   /* K = 1.0 / K; */
   /* omega = 1.0 / (omega-0.5); */
   /* omega = (omega-0.5) / K; */
   /* omega = (omega-0.5) / sqrt(K);
   omega += 0.5; */

   omega = sqrt (omega / K);
   omega += 0.5;

   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      x = K * (x-omega) * (1.0 - omega - x); 
      x -= floor (x);
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a 
       running sum). */
   xpoint = x;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      x = K * (x-omega) * (1.0 - omega - x); 
      x -= floor (x);

      y = fabs (x-xpoint);
      y -= floor (y);
      if (y < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (NORMLOGIS * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the folded map
 */

float folded_poincare_recurrance_time (CircleData *dat,
                       double omega,
                       double K)

{
   double	x, y;
   double	xpoint;
   int		iter;
   long		num_recurs, time_recur=0;
   
   /* First, we give a spin for 500 cycles, giving the non-chaotic 
    * parts a chance to phase-lock */
   x = 0.0;
   for (iter=0; iter<SETTLE_TIME; iter++) {
      x -= omega;
      if (x > 0.5) {
         x = 1.0 - x;
      } 
      x *= K;
      x -= floor (x);
   }

   /* OK, now, we begin to measure the average amount of time to recur */
   /* (note that we don't have todo += with iter, since its already a running sum). */
   xpoint = x;
   num_recurs = 0;
   for (iter=0; iter < dat->itermax; iter++) {
      x -= omega;

      /*   x = K * x * (1.0 - x);  -- this is logistic */
      if (x > 0.5) {
         x = 1.0 - x;
      } 
      x *= K;
      x -= floor (x);

      y = fabs (x-xpoint);
      y -= floor (y);
      if (y < EPSILON) {
         num_recurs ++;
         time_recur = iter;
      }
   }

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / (NORMALIZATION * ((double)num_recurs));

   return ((float) x);
}

/*-------------------------------------------------------------------*/
   

void do_circle (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;

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
   omega_max = 1.0;
   /* K_min = 1.0; */
   /* special for irs printer */
   K_min = 2500.0 / 1600.0;
   K_max = 0.0;

   /* Close up of central flame */
   omega_min = 5.0/12.0;
   omega_max = 7.0/12.0;
   K_min = 1.0/3.0;
   K_max = 1.0/6.0;

   omega_min = 0.5 - 0.05;
   omega_max = 0.5 + 0.05;
   K_min = 1.0/3.0 + 0.05;
   K_max = 1.0/3.0 - 0.05;

   omega_min = -0.5;
   omega_max = 0.5;
   K_min = 1.2;
   K_max = 0.2;

   omega_min = 0.25;
   omega_max = 0.75;
   K_min = 1.25;
   K_max = 0.75;

   omega_min = 0.0;
   omega_max = 1.0;
   K_min = 3.00;
   K_max = 0.0;

   dat.itermax = 2000;
   
   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              (DensityCB) circle_poincare_recurrance_time, &dat);

#ifdef TONGUES_TRI
   /* closeup of Sinai's Toungues */
   omega_min = 0.0;
   omega_max = 1.0;
   K_min = 0.5;
   K_max = 0.0;
   walk_utri (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              circle_poincare_recurrance_time, &dat);
#endif

   
   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "standard circle map\n");
   fprintf (fp, "x(n+1) = x(n) + Omega - K * sin(2*PI*x(n))\n");
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", NORMALIZATION);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "K on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);
   
}

/*-------------------------------------------------------------------*/
   
void do_im_circle (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   /*   The BIG picture
   */
   omega_min = 0.1;
   omega_max = 0.106;
   /* K_min = 1.0; */
   /* special for irs printer */
   K_min = 1.5;
   K_max = 0.0;

   dat.itermax = 11000;
   
   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              (DensityCB) im_circle_poincare_recurrance_time, &dat);

   
   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "imaginary K coupling circle map\n");
   fprintf (fp, "x(n+1) = x(n) + Omega - K * sin(2*PI*x(n))\n");
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", INOR);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "K on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);
   
}

/*-------------------------------------------------------------------*/
   
void do_co_circle (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		iii, itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;
   char buff[100];

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
   omega_max = 1.0;
   /* K_min = 1.0; */
   /* special for irs printer */
   K_min = 3.0;
   K_max = 0.0;

   dat.itermax = 6000;
   
/*
complex_phase =  -0.0001;
for (iii=0; iii<11; iii++) {
complex_phase +=0.0001;
*/
{
iii=4;
complex_phase = 0.0004;

   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              (DensityCB) co_circle_poincare_recurrance_time, &dat);

#ifdef TONGUES_TRI
   /* closeup of Sinai's Toungues */
   omega_min = 0.0;
   omega_max = 1.0;
   K_min = 0.5;
   K_max = 0.0;
   walk_utri (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              circle_poincare_recurrance_time, &dat);
#endif

   
   /* dump the floating point data */
   sprintf (buff, "%s%d", filename, iii);
   if ( (fp = Fopen (buff, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (buff, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "Complex K coupling circle map\n");
   fprintf (fp, "Phase of K is %f \n", complex_phase);
   fprintf (fp, "x(n+1) = x(n) + Omega - K * sin(2*PI*x(n))\n");
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", CNOR);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "K on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);

}
   
}

/*-------------------------------------------------------------------*/
   
void do_logistic (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   /*   The BIG picture
   */
   omega_min = 0.66666;
   omega_max = 1.66666;
   K_min = 10.0;
   K_max = 0.0;
   
   dat.itermax = 15000;

   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              /* clamp_logistic_poincare_recurrance_time, &dat); */
              (DensityCB) squash_logistic_poincare_recurrance_time, &dat);

   
   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "standard logistic map\n");
   fprintf (fp, "y = x(n) - (sqrt(omega/K)+0.5) \n");
   fprintf (fp, "z = x(n) + (sqrt(omega/K)+0.5) \n");
   /* fprintf (fp, "y = x(n) - (omega-0.5)/sqrt(K) - 0.5\n");  */
   fprintf (fp, "z = x(n) + (omega-0.5)/sqrt(K) + 0.5\n"); 
   /* fprintf (fp, "y = x(n) - Omega \n");  */
   fprintf (fp, "z = x(n) + Omega \n"); 
   /* fprintf (fp, "x(n+1) = lambda * y * (1-z) mod 1 \n"); */
   /* fprintf (fp, "x(n+1) = (1/lambda) * y * (1-z) mod 1 \n"); */
   /* fprintf (fp, "x(n+1) = lambda * y * (1-y) mod 1 \n"); */
   /* fprintf (fp, "x(n+1) = (lambda * sin (pi*y)) mod 1 \n"); */
   /* fprintf (fp, "x(n+1) = (lambda * sin (pi*y) / (pi*y)) mod 1 \n"); */
   /* fprintf (fp, "x(n+1) = (lambda * exp (y*y)) mod 1 \n"); */
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", NORMLOGIS);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "lambda on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);
   
}

/*-------------------------------------------------------------------*/
   
void do_folded (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		itermax=0;
   float	omega_min, omega_max, K_min, K_max;
   FILE		*fp;
   CircleData dat;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   /*   The BIG picture
   */
   omega_min = -1.0;
   omega_max = 1.0;
   /* K_min = 14.0; */
   K_min = 4.0; 
   K_max = 0.0;
   
   dat.itermax = 5000;

   /* fill it in */
   walk_rect (data, data_width, data_height,
              omega_min, omega_max, K_min, K_max, 
              (DensityCB) folded_poincare_recurrance_time, &dat);

   
   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".text")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "standard folded map\n");
   fprintf (fp, "y = x(n) - Omega \n");
   fprintf (fp, "x(n+1) = lambda * (0.5 - abs(y-0.5)) mod 1 \n");
   fprintf (fp, "Each pixel shows number of cycles to recur within EPSILON \n");
   fprintf (fp, "where EPSILON = %f \n", EPSILON);
   fprintf (fp, "But first, system allowed to settle for %d cycles\n", SETTLE_TIME);
   fprintf (fp, "Recurrance time normalized by %f\n", NORMALIZATION);
   fprintf (fp, "\n");
   fprintf (fp, "floating point values, one per pixel, (0<x<1); \n");
   fprintf (fp, "Omega on X axis min = %f max = %f \n", omega_min, omega_max);
   fprintf (fp, "lambda on Y axis min = %f max = %f \n", K_max, K_min);
   fprintf (fp, "width: %d height: %d \n", data_width, data_height);
   fprintf (fp, "Iteration stopped after %d iterations \n", itermax);
   fclose (fp);
   
}

/*-------------------------------------------------------------------*/
#include <signal.h>
#include <unistd.h>

int 
main (int argc, char *argv[]) 

{

   if (4 > argc) {
      fprintf (stderr, "Usage: %s <out-file> <width> <height> [<param>] \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   signal (SIGFPE, SIG_IGN);
   if (argc >4) param = atoi (argv[4]);

   if (!strcmp ("moire", argv[0])) do_moire (argv[1], atoi (argv[2]), atoi (argv[3]));

   /* do_circle (argv[1], atoi (argv[2]), atoi (argv[3])); */
   /* do_im_circle (argv[1], atoi (argv[2]), atoi (argv[3])); */
   /* do_co_circle (argv[1], atoi (argv[2]), atoi (argv[3])); */
   /* do_logistic (argv[1], atoi (argv[2]), atoi (argv[3])); */
   /* do_folded (argv[1], atoi (argv[2]), atoi (argv[3])); */

   if (!strcmp ("mandel", argv[0])) do_mandel (argv[1], atoi (argv[2]), atoi (argv[3])); 
   if (strstr (argv[0], "tong")) do_cf_gap (argv[1], atoi (argv[2]), atoi (argv[3])); 

   return 0;
}

/*-------------------- END OF FILE -----------------------------------*/
