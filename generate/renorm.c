/*
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "opers.h"

/*-------------------------------------------------------------------*/
extern FILE *Fopen();
extern FILE *Fopenr();

void main (int argc, char *argv[]) 
{
   float	*data_in ;	/* my data array */
   // float	*data_out;	/* my data array */
   // float		*tmp; 		/* pointer into data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		i,j;
   char 		str[80];
   FILE		*fp_in, *fp_out;
   float 		scale_fact;
   float 		offset;
   
   if (argc < 3) {
      fprintf (stderr, "usage: %s <input file> <output file> [<scale factor>] [<offset>] \n", argv[0]);
      return;
   }

   fprintf (stderr, "%s %s %s ", argv[0], argv[1], argv[2]);
   
   if (4 <= argc) {
      scale_fact = (float) atof (argv[3]);
      fprintf (stderr, "%g ", scale_fact);
   } else {
      scale_fact = 1.0;
   }

   if (5 <= argc) {
      offset = (float) atof (argv[4]);
      fprintf (stderr, "%g ", offset);
   } else {
      offset = 0.0;
   }
   fprintf (stderr, "\n");

   
   /*-----------------------------------------------*/
   /* open input file */
   if ( (fp_in = Fopenr (argv[1], ".flo")) == NULL) {
      fprintf (stderr, " Can't open input file %s \n", argv[1]);
      return;
   }

   i = 0;
   j = 1;
   while ((0 != j) && ('\n' != j)) {
      j =  fgetc (fp_in);
      str[i] = j;
      i++;
   }
   
/*
   if ( NULL == fgets (str, 80, fp_in)) {
      fprintf (stderr, " Can't read input file %s \n", argv[1]);
      fclose (fp_in);
      return;
   }
*/
   
   sscanf (str, "%d %d", &data_width, &data_height);
   fprintf (stderr, "Input file %s has width %d height %d \n", argv[1], data_width, data_height);
   
   data_in = (float *) malloc (data_width * data_height * sizeof(float));
   
   /* read floating point data */
   fread(data_in, sizeof(float), data_width*data_height, fp_in);
   fclose (fp_in);
   
   /*-----------------------------------------------*/
   /* open output file */
   if ( (fp_out = Fopen (argv[2], ".flo")) == NULL) {
      fprintf (stderr, " Can't open output file %s \n", argv[2]);
      return;
   }
   
   /* dump the size to output */
   fprintf (fp_out, "%d %d\n", data_width, data_height);
   
   /*-----------------------------------------------*/
   
   if (!strcmp ("takelog", argv[0])) takelog (data_in, data_width, data_height);
   if (!strcmp ("abs", argv[0])) absolval (data_in, data_width, data_height);
   if (!strcmp ("renorm", argv[0])) expand (data_in, data_width, data_height, scale_fact, offset);
   if (!strcmp ("reclamp", argv[0])) rescale (data_in, data_width, data_height, scale_fact);
   if (!strcmp ("clamp", argv[0])) clamp (data_in, data_width, data_height, scale_fact, offset);
   if (!strcmp ("dump", argv[0])) dump (data_in, data_width, data_height, scale_fact);

/*
   aver = avg (data_in, data_width, data_height);
   msq = sqdev (data_in, data_width, data_height);
   msq = sqrt (msq);
   printf ("yo average was %f rms %f \n", aver, msq);
   expand (data_in, data_width, data_height, scale_fact/msq, aver);
*/
   
   fwrite (data_in, sizeof(float), data_width*data_height, fp_out);
   
   fclose (fp_out);
}

/* --------------- END OF FILE ------------------------- */
