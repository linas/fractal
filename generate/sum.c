/*
 * FUNCTION:
 * take sum of two things
 *
 * HISTORY:
 * Created Linas Vepstas October 1989
 * Updates March 1996
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "opers.h"

/*-------------------------------------------------------------------*/

void sum (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      a[i] += b[i];
   }
}

/*-------------------------------------------------------------------*/

void diff (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      a[i] -= b[i];
   }
}

/*-------------------------------------------------------------------*/

void paste (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      if (b [i] > 0.01) a[i] = b[i];
   }
}

/*-------------------------------------------------------------------*/

void angle (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      a[i] = M_PI + atan2 (a[i], b[i]);
      a[i] *= 0.5 / M_PI;
   }
}

/*-------------------------------------------------------------------*/

void curl (float a[], float b[], int nx, int ny) {
   int i, j;
   for (j=0; j<ny-1; j++) {
      for (i=0; i<nx-1; i++) {
         a[nx*j+i] -= a[nx*(j+1)+i];
         b[nx*j+i] -= b[nx*j+i+1];
         a[nx*j+i] = b[nx*j+i] - a[nx*j+i];
      }
   }
}

/*-------------------------------------------------------------------*/

void divergence (float a[], float b[], int nx, int ny) {
   int i, j;
   for (j=0; j<ny-1; j++) {
      for (i=0; i<nx-1; i++) {
         a[nx*j+i] -= a[nx*j+i+1];
         b[nx*j+i] -= b[nx*(j+1)+i];
         a[nx*j+i] += b[nx*j+i];
      }
   }
}

/*-------------------------------------------------------------------*/

extern FILE *Fopen();
extern FILE *Fopenr();

int 
main (int argc, char *argv[]) 

{
   float	*data_a;		/* my data array */
   float	*data_b;		/* my data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		globlen;
   char 	str[80];
   FILE		*fp_a, *fp_b, *fp_out;
   
   if (argc < 4) {
      printf ("Usage: %s <input file> <input file> <output file>\n", argv[0]);
      return 1;
   }

   /*----------------------------------------------------*/
   /* open input file */
   if ( (fp_a = Fopenr (argv[1], ".flo")) == NULL) {
      printf (" Can't open input file %s \n", argv[1]);
      return 1;
   }
   if ( NULL == fgets (str, 80, fp_a)) {
      printf (" Can't read input file %s \n", argv[1]);
      fclose (fp_a);
      return 1;
   }

   sscanf (str, "%d %d", &data_width, &data_height);
   printf (" input file %s has width %d height %d \n", argv[1], data_width, data_height);
   
   /*----------------------------------------------------*/
   /* open input file */
   if ( (fp_b = Fopenr (argv[2], ".flo")) == NULL) {
      printf (" Can't open input file %s \n", argv[2]);
      return 1;
   }
   if ( NULL == fgets (str, 80, fp_b)) {
      printf (" Can't read input file %s \n", argv[2]);
      fclose (fp_b);
      return 1;
   }

   sscanf (str, "%d %d", &data_width, &data_height);
   printf (" input file %s has width %d height %d \n", argv[2], data_width, data_height);
   
   /*----------------------------------------------------*/
   
   globlen = data_width*data_height;
   data_a = (float *) malloc (globlen * sizeof (float));
   data_b = (float *) malloc (globlen * sizeof (float));
   
   /* read floating point data */
   fread(data_a, sizeof(float), globlen, fp_a);
   fread(data_b, sizeof(float), globlen, fp_b);
   fclose (fp_a);
   fclose (fp_b);

   /*----------------------------------------------------*/
   if (!strcmp (argv[0], "sum")) sum (data_a, data_b, globlen);
   if (!strcmp (argv[0], "diff")) diff (data_a, data_b, globlen);
   if (!strcmp (argv[0], "paste")) paste (data_a, data_b, globlen);
   if (!strcmp (argv[0], "angle")) angle (data_a, data_b, globlen);
   if (!strcmp (argv[0], "curl")) curl (data_a, data_b, data_width, data_height);
   if (!strcmp (argv[0], "div")) divergence (data_a, data_b, data_width, data_height);

   /*----------------------------------------------------*/
   /* open output file */
   if ( (fp_out = Fopen (argv[3], ".flo")) == NULL) {
      printf (" Can't open output file \n");
      return 1;
   }

   /* dump the size to output */
   fprintf (fp_out, "%d %d\n", data_width, data_height);
   fwrite (data_a, sizeof(float), globlen, fp_out);
   fclose (fp_out);

   free (data_a);
   free (data_b);
   
   return 0;
}

/* --------------------------- END OF FILE ------------------------- */
