/*
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

#include <stdio.h>
#include <math.h>
#include "opers.h"

/*-------------------------------------------------------------------*/
extern FILE *Fopen();
extern FILE *Fopenr();

void main (int argc, char *argv[]) 
{
   float	*data_in ;	/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		i,j;
   char 	str[80];
   FILE		*fp_in; 
   double 	aver, msq;
   
   if ( argc < 2) {
      printf ("usage: %s <input file> \n", argv[0]);
      return;
      }
   
   /*-----------------------------------------------*/
   /* open input file */
   if ( (fp_in = Fopenr (argv[1], ".flo")) == NULL) {
      printf (" Can't open input file \n");
      return;
      }
   if ( NULL == fgets (str, 80, fp_in)) {
      printf (" Can't read input file \n");
      fclose (fp_in);
      return;
      }
   
   sscanf (str, "%d %d", &data_width, &data_height);
   printf (" input file has width %d height %d \n", data_width, data_height);
   
   data_in = (float *) malloc (data_width * data_height * sizeof(float));
   
   /* read floating point data */
   fread(data_in, sizeof(float), data_width*data_height, fp_in);
   fclose (fp_in);
   
   aver = avg (data_in, data_width, data_height);
   msq = sqdev (data_in, data_width, data_height);
   msq = sqrt (msq);
   printf ("yo average was %f rms %f \n", aver, msq);
}

/* --------------- END OF FILE ------------------------- */
