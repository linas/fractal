/*
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "opers.h"
#include "fileio.h"

/*-------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
   float 		scale_fact;
   float 		offset;

   if (argc < 3) {
      fprintf (stderr, "Usage: %s <input file> <output file> [<scale factor>] [<offset>] \n", argv[0]);
      return 1;
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
	const char* suff;
   unsigned int data_width, data_height; /* data array dimensions */
	float* data_in = read_floats(argv[1], &suff, &data_width, &data_height);
	if (NULL == data_in)
	{
		fprintf (stderr, " Can't open input file %s \n", argv[1]);
		return 2;
	}

   /*-----------------------------------------------*/
	/* Strip out the path */
	char * p = strdup(argv[0]);
	char * name = strrchr(p, '/');
	if (NULL == name) name = p;
	else name++;

   if (!strcmp ("takelog", name)) takelog (data_in, data_width, data_height);
   else if (!strcmp ("recip", name)) recip (data_in, data_width, data_height);
   else if (!strcmp ("abs", name)) absolval (data_in, data_width, data_height);
   else  if (!strcmp ("takeroot", name)) takeroot (data_in, data_width, data_height, scale_fact);
   else  if (!strcmp ("renorm", name)) expand (data_in, data_width, data_height, scale_fact, offset);
   else  if (!strcmp ("reclamp", name)) rescale (data_in, data_width, data_height, scale_fact);
   else  if (!strcmp ("clamp", name)) clamp (data_in, data_width, data_height, scale_fact, offset);
   else  if (!strcmp ("dump", name)) dump (data_in, data_width, data_height);
	else
	{
		fprintf(stderr, "Error: unknown transform %s\n", name);
		exit(1);
	}

/*
   aver = avg (data_in, data_width, data_height);
   msq = sqdev (data_in, data_width, data_height);
   msq = sqrt (msq);
   printf ("yo average was %f rms %f \n", aver, msq);
   expand (data_in, data_width, data_height, scale_fact/msq, aver);
*/

	write_floats(argv[2], suff, data_in, data_width, data_height);

	return 0;
}

/* --------------- END OF FILE ------------------------- */
