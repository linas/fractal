/* 
 * NAME:
 * gap-tongue.C
 *
 * FUNCTION:
 * generate pixmaps showing continued fraction gaps
 *
 * HISTORY:
 * created Linas Vepstas Jan 1991
 * spruced up Linas Vepstas July 1993
 * fixed major performance bug (AIX fmod() stinks!) -- January 1994
 * added logistic map -- February 1994
 * added classic mandelbrot -- June 1995
 * Added continued fraction gaps -- October 2004
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Farey.h"
#include "image.h"
#include "util.h"

/*-------------------------------------------------------------------*/

typedef struct gap_data_s {
	ContinuedFraction f;
} GapData;

/*-------------------------------------------------------------------*/
/*
 * Compute the wavy gaps for a given continued fraction.
 */

float cf_gap (GapData *dat,
                       double x,
                       double w)

{
	dat->f.SetReal (x);
	double g = dat->f.ToXPlus (w);
	return g;
}

/*-------------------------------------------------------------------*/
   
void do_cf_gap (char * filename, int width, int height)

{
   float	*data ;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   float	x_min, x_max, w_min, w_max;
   FILE		*fp;
   GapData dat;

   data = (float *) malloc (sizeof (float) * width * height);
   if (data == NULL) {
      fprintf (stderr, "Whoops, no memory \n");
      exit (EXIT_FAILURE);
   }

   data_width = width;
   data_height = height;
   /*   The BIG picture
   */
   x_min = 0.0;
   x_max = 1.0;
   w_min = 0.0; 
   w_max = 1.0;

   
   /* fill it in */
   walk_rect (data, data_width, data_height,
              x_min, x_max, w_min, w_max, 
              (DensityCB) cf_gap, &dat);

   
   /* dump the floating point data */
   if ( (fp = Fopen (filename, ".flo")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fp);
   fclose (fp);
   
   if ( (fp = Fopen (filename, ".txt")) == NULL) {
      fprintf (stderr, " File open failure \n");
      exit (EXIT_FAILURE);
   }

   fprintf (fp, "gappy stuff\n");
   fclose (fp);
}

/* ----------------------- end of file  ---------------------- */
   
