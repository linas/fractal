
/* 
 * HISTORY:
 * Quick and dirty hack -- linas November 1989 
 * Quick and dirty updates -- linas June 1991 
 * Refurbish, dust off, -- linas March 1996 
 */

#include <stdio.h>
#include <math.h>

#include "opers.h"


/*-------------------------------------------------------------------*/

void edge (row_m1, row, row_p1, row_out, sizex, scale_factor)
float  		row_m1 []; /* row minus 1 */
float  		row [];    /* row */
float  		row_p1 []; /* row plus 1 */
float  		row_out [];
unsigned int	sizex;
float		scale_factor;
{
   long		i;
   float 		oops;
   
   oops = 1.0 / scale_factor;
   /* first pixel is special */
   row_out[0] = 3.0*row[0] -row[1] -row_p1[0] -row_m1[0];
   row_out[0] *= oops;
   
   for (i=1; i<sizex-1; i++) {
      row_out[i] = 4.0*row[i] -row[i+1] -row[i-1] -row_p1[i] -row_m1[i];
      row_out[i] *= oops;
   }

   /* last pixel is special */
   row_out[sizex-1] = 3.0*row[sizex-1] -row[sizex-2] -row_p1[sizex-1] -row_m1[sizex-1];
   row_out[sizex-1] *= oops;
}

/*-------------------------------------------------------------------*/

void dredge (row_m1, row, row_p1, row_out, sizex, threshold)
float  		row_m1 []; /* row minus 1 */
float  		row [];    /* row */
float  		row_p1 []; /* row plus 1 */
float  		row_out [];
unsigned int	sizex;
float		threshold;
{
   long		i;
   
   /* first pixel is special */
   if ((fabs (row[0]-row[1]) > threshold)  ||
       (fabs (row_p1[0]-row[0]) > threshold)  ||
       (fabs (row_m1[0]-row[0]) > threshold)) {
       row_out [0] = 0.999;
    } else {
       row_out[0] = 0.0;
    }

   for (i=1; i<sizex-1; i++) {
      if ((fabs (row[i]-row[i-1]) > threshold)  ||
          (fabs (row[i+1]-row[i]) > threshold)  ||
          (fabs (row_p1[i]-row[i]) > threshold)  ||
          (fabs (row_m1[i]-row[i]) > threshold)) {
           row_out[i] = 0.999;
        } else {
           row_out[i] = 0.0;
        }
   }

   /* last pixel is special */
   if ((fabs (row[sizex-1]-row[sizex-2]) > threshold)  ||
       (fabs (row_p1[sizex-1]-row[sizex-1]) > threshold)  ||
       (fabs (row_m1[sizex-1]-row[sizex-1]) > threshold)) {
       row_out[sizex-1] = 0.999;
    } else {
       row_out[sizex-1] = 0.0;
    }
}

/*-------------------------------------------------------------------*/

void contour (row_m1, row, row_p1, row_out, sizex)
float  		row_m1 []; /* row minus 1 */
float  		row [];    /* row */
float  		row_p1 []; /* row plus 1 */
float  		row_out [];
unsigned int	sizex;
{
   long		i;
   
   /* first pixel is special */
   if ( (((int) row[0]) != ((int) row[1]))  ||
        (((int) row[0]) != ((int) row_p1[0]))  ||
        (((int) row[0]) != ((int) row_m1[0])) ) {
        row_out[0] = 0.999;
     } else {
        row_out[0] = 0.0001;
     }

   for (i=1; i<sizex-1; i++) {
      if ( (((int) row[i]) != ((int) row[i-1]))  ||
           (((int) row[i]) != ((int) row[i+1]))  ||
           (((int) row[i]) != ((int) row_p1[i]))  ||
           (((int) row[i]) != ((int) row_m1[i])) ) {
           row_out[i] = 0.999;
        } else {
           row_out[i] = 0.0001;
        }
   }

   /* last pixel is special */
   if ( (((int) row[sizex-1]) != ((int) row[sizex-2]))  ||
        (((int) row[sizex-1]) != ((int) row_p1[sizex-1]))  ||
        (((int) row[sizex-1]) != ((int) row_m1[sizex-1])) ) {
        row_out[sizex-1] = 0.999;
     } else {
        row_out[sizex-1] = 0.0001;
     }

}

/*-------------------------------------------------------------------*/

#define SCALE 1.3
#define THRESH 0.9

#define PIX_WID 1600
#define PIX_HIG 1280
extern FILE *Fopen();
extern FILE *Fopenr();

void main (argc, argv) 
int argc;
char *argv[];

{
   static float	data_in [PIX_WID*4];	/* my data array */
   static float	data_out [PIX_WID*4];	/* my data array */
   float		*ptr_a; 		/* pointer into data array */
   float		*ptr_b; 		/* pointer into data array */
   float		*ptr_c; 		/* pointer into data array */
   float		*tmp; 		/* pointer into data array */
   unsigned int	data_width, data_height;/* data array dimensions */
   int		i,j;
   char 		str[80];
   FILE		*fp_in, *fp_out;
   float 		scale_fact;
   
   if ( argc < 3) {
      printf ("usage: %s <input file> <output file> \n", argv[0]);
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
   
   /*-----------------------------------------------*/
   /* open output file */
   if ( (fp_out = Fopen (argv[2], ".flo")) == NULL) {
      printf (" Can't open output file \n");
      return;
   }

   /* dump the size to output */
   fprintf (fp_out, "%d %d\n", data_width, data_height);
   
   /*-----------------------------------------------*/
   ptr_a = data_in;
   ptr_b = ptr_a + data_width;
   ptr_c = ptr_b + data_width;
   /* read in first two rows */
   fread(ptr_a, sizeof(float), data_width, fp_in);
   fread(ptr_b, sizeof(float), data_width, fp_in);
   
   /* first row is special */
   expand (ptr_a, data_width, 1, SCALE, 0.0);
   expand (ptr_b, data_width, 1, SCALE, 0.0);
   /*
   dredge (ptr_a, ptr_b, ptr_b, data_out, data_width, THRESH);
   */
   contour (ptr_a, ptr_b, ptr_b, data_out, data_width);

   fwrite (data_out, sizeof(float), data_width, fp_out);
   
   for (i=1; i<data_height-1; i++) {
       printf (" working on row %d\n", i);
   
       /* read floating point data */
       fread(ptr_c, sizeof(float), data_width, fp_in);

       expand (ptr_c, data_width, 1, SCALE, 0.0);

       contour (ptr_a, ptr_b, ptr_c, data_out, data_width);
/*
       dredge (ptr_a, ptr_b, ptr_c, data_out, data_width, THRESH);
*/
   
       /* dump the floating point data */
       fwrite (data_out, sizeof(float), data_width, fp_out);
   
       /* do the Harlem shuffle */
       tmp = ptr_a;
       ptr_a = ptr_b;
       ptr_b = ptr_c;
       ptr_c = tmp;
    }

   /* last row is special */
   /*
   dredge (ptr_b, ptr_b, ptr_c, data_out, data_width, THRESH);
   */
   contour (ptr_b, ptr_b, ptr_c, data_out, data_width);
   fwrite (data_out, sizeof(float), data_width, fp_out);
   
   fclose (fp_in);
   fclose (fp_out);

}

/* --------------- END OF FILE ------------------------- */
