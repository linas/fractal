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
#include <string.h>

#include "opers.h"
#include "fileio.h"

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

void mask (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      if (b[i] > 0.5) a[i] = 0.0;
   }
}

/*-------------------------------------------------------------------*/

void paste (float a[], float b[], int len) {
   int i;
   for (i=0; i<len; i++) {
      if (b[i] > 0.01) a[i] = b[i];
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

int main (int argc, char *argv[])

{
   if (argc < 4) {
      printf ("Usage: %s <input file> <input file> <output file>\n", argv[0]);
      return 1;
   }

   /*----------------------------------------------------*/
   /* open input file */
	const char* suff;
   unsigned int	data_width, data_height;/* data array dimensions */
   float* data_a = read_floats(argv[1], &suff, &data_width, &data_height);
   printf ("input file %s.%s has width %d height %d\n",
          argv[1], suff, data_width, data_height);

   float* data_b = read_floats(argv[2], &suff, &data_width, &data_height);
   printf ("input file %s.%s has width %d height %d\n",
          argv[2], suff, data_width, data_height);

	int globlen = data_width*data_height;

   /*----------------------------------------------------*/
   /* strip out leading pathame */
   char	*cmd = rindex (argv[0], '/');
   if (!cmd) {
      cmd = argv[0];
   } else {
      cmd++;
   }

   if (!strcmp (cmd, "sum")) sum (data_a, data_b, globlen);
   if (!strcmp (cmd, "diff")) diff (data_a, data_b, globlen);
   if (!strcmp (cmd, "mask")) mask (data_a, data_b, globlen);
   if (!strcmp (cmd, "paste")) paste (data_a, data_b, globlen);
   if (!strcmp (cmd, "angle")) angle (data_a, data_b, globlen);
   if (!strcmp (cmd, "curl")) curl (data_a, data_b, data_width, data_height);
   if (!strcmp (cmd, "div")) divergence (data_a, data_b, data_width, data_height);

   /*----------------------------------------------------*/
   /* open output file */
	write_floats(argv[3], suff, data_a, data_width, data_height);

   return 0;
}

/* --------------------------- END OF FILE ------------------------- */
