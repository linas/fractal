/* 
 * NAME:
 * draw.c
 *
 * FUNCTION:
 * generate pixmaps of circle
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

#include "Farey.h"

extern FILE *Fopen();

/*-------------------------------------------------------------------*/
/* initialize drawing area */

void 
initialize (float *glob, unsigned int sizex, unsigned int sizey)
{
   unsigned int	i, globlen;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
}

/*-------------------------------------------------------------------*/
/* draws a circle */

void draw_circle (float *glob, unsigned int sizex, unsigned int sizey,
                  double centerx, double centery, double radius)
{
   unsigned int	i;
   int		ic, is;
   double theta, dtheta;
   double tmp, ds, dc, s, c;
   
   dtheta = 1.0 / radius;
   ds = sin (dtheta);
   dc = cos (dtheta);
   s = radius;
   c = 0.0;
   theta = 0.0;

   for (i=0; i< 7.0*radius; i++) {
      ic = (int) (c + centerx);
      is = (int) (s + centery);
      if (0 <= ic  && ic < sizex && 0 <= is && is < sizey) {
         glob [is*sizex + ic] = 1.0;
      }
      tmp = s*dc + c*ds;
      c = c*dc - s*ds;
      s = tmp;
      theta += dtheta;
   }
}

/*-------------------------------------------------------------------*/
/* draws a curve in polar coordinates */
typedef double (*Param)(double);

void draw_parameter (float *glob, unsigned int sizex, unsigned int sizey,
                  double centerx, double centery, double scale, 
                  Param f, Param g)
{
   unsigned int	i,imax;
   double t, r, theta;
   int		ix, iy;
   
   imax = sizex *10;
   for (i=0; i< imax; i++) {
      t = ((double) i) / ((double) imax);

      theta = f(t);
      r = scale * g(t);
      
      ix = (int) (centerx + r * cos (theta));
      iy = (int) (centery + r * sin (theta));
      if (0 <= ix  && ix < sizex && 0 <= iy && iy < sizey) {
         glob [iy*sizex + ix] = 1.0;
      }
   }
}

/*-------------------------------------------------------------------*/

double thet (double t) { return 2.0*M_PI*t; }
double rad (double t) {return sin (M_PI*t); }

struct Farey *far = NULL;

/* a Farey distortion */
double fthet (double t) 
{ 
   double ecc;

   /* initialize Farey number class */
   if (!far) far = CreateFarey();

   SetReal (far, t);
   ecc = ContinuedFractionToFarey (far);
   return 2.0*M_PI*ecc; 
}


double frad (double t) 
{ 
   double ecc;

   /* initialize Farey number class */
   if (!far) far = CreateFarey();

   SetReal (far, t);
   ecc = ContinuedFractionToFarey (far);
   return sin (M_PI * ecc); 
}


/*-------------------------------------------------------------------*/
/* draws the mandelbrot bulb */

void 
draw_bulb (float *glob, unsigned int sizex, unsigned int sizey,
                  double centerx, double centery, double scale)
{
   unsigned int	i,imax;
   double t, r, theta;
   int		ix, iy;
   double ecc;

   /* initialize Farey number class */
   if (!far) far = CreateFarey();
   
   imax = sizex *10;
   for (i=0; i< imax; i++) {
      t = ((double) i) / ((double) imax);

      /* cardiod in polar coords */
      theta = 2.0*M_PI*t;
      r = scale *  sin (M_PI*t);

      /* bud size parameter */
      SetReal (far, t);
      ecc = ContinuedFractionToFarey (far);

      /* lay down a point */
      ix = (int) (centerx + r * cos (theta));
      iy = (int) (centery + r * sin (theta));
      if (0 <= ix  && ix < sizex && 0 <= iy && iy < sizey) {
         glob [iy*sizex + ix] = 1.0;
      }

      draw_circle (glob, sizex, sizey, ix, iy, scale*ecc);
   }
}

/*-------------------------------------------------------------------*/
/* draws a median */

void draw_median_circle (float *glob, unsigned int sizex, unsigned int sizey,
                  double centerx, double centery, double radius)
{
   unsigned int	i, imax;
   int		ic, is;
   double theta, dtheta;
   double rm, a,g, s, c;
   
   s = radius;
   c = 0.0;

   imax = 7.0*radius;
   for (i=0; i<imax; i++) {
      theta =  2.0 * M_PI * ((double) i) / ((double) imax);

      c = 1.01 + cos (theta);
      s = 1.01 + sin (theta);

      a = 0.5 * (c+s);
      g = sqrt (c*s);
      rm = s/c;

      a *= radius;
      g *= radius;

      ic = (int) (a + centerx);
      is = (int) (g + centery);
      if (0 <= ic  && ic < sizex && 0 <= is && is < sizey) {
         glob [is*sizex + ic] = 1.0;
      }
   }
}

/*-------------------------------------------------------------------*/
   
void do_circle (char * filename, int width, int height,
                double cx, double cy, double scale)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   FILE		*fp;
   char buff[100];

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   initialize (data, data_width, data_height);

   // draw_parameter (data, data_width, data_height, cx, cy, scale, thet, rad);
   // draw_bulb (data, data_width, data_height, cx, cy, scale);
   draw_median_circle (data, data_width, data_height, cx, cy, scale);
   
   /* dump the floating point data */
   sprintf (buff, "%s", filename);
   if ( (fp = Fopen (buff, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
}

/*-------------------------------------------------------------------*/

int 
main (int argc, char *argv[]) 

{

   if (argc != 7) {
      fprintf (stderr, "Usage: %s <out-file> <width> <height> <centerx> <centery> <radius>\n", argv[0]);
      exit (EXIT_FAILURE);
   }

   do_circle (argv[1], atoi (argv[2]), atoi (argv[3]),
              atof(argv[4]), atof(argv[5]), atof(argv[6])); 

   exit (EXIT_SUCCESS);
}

/*-------------------- END OF FILE -----------------------------------*/
