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

/*-------------------------------------------------------------------*/
/* draws a circle */

void draw_circle (float *glob, unsigned int sizex, unsigned int sizey)
{
   unsigned int	i,j, globlen;
   double x, tpx, tpy, phi, car, cir, card;
   int icar, icir, icard, jcard;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

   for (i=0; i< sizex; i++) {

      x = (double) i / (double) sizex;

      cir = 2.0*x -1;
      cir = -0.5 + sqrt (25.0/16.0 - cir*cir);

      car = 0.5 * sqrt (1.25 - cos (2.0*M_PI*x));

      x *= 2.0 * M_PI;
      tpy = sin (x) - 0.5 * sin(2.0*x);
      tpx = cos (x) - 0.5 * cos(2.0*x);
      phi = atan2 (tpy, tpx);
      card =  0.5 * sqrt (1.25 - cos (phi));

      phi /= 2.0 * M_PI;
      jcard = sizex * phi;
      if (0 > jcard) jcard = - jcard;

printf ("its %d %d i, jcard\n", i, jcard);

      icar = sizex - sizex*0.5*car;
      icard = sizex - sizex*0.5*card;
      icir = sizex - sizex*0.5*cir;

      glob [icir*sizex + i] = 0.5;
      glob [icar*sizex + jcard] = 0.8;
      glob [icar*sizex + i] = 1.0;
   }
}

/*-------------------------------------------------------------------*/
   
void do_circle (char * filename, int width, int height)

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

   draw_circle (data, data_width, data_height);
   
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

void main (int argc, char *argv[]) 

{

   if (argc != 4) {
      fprintf (stderr, "Usage: %s <out-file> <width> <height>\n", argv[0]);
      exit (EXIT_FAILURE);
   }

   do_circle (argv[1], atoi (argv[2]), atoi (argv[3]));

   exit (EXIT_SUCCESS);
}

/*-------------------- END OF FILE -----------------------------------*/
