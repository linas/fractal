
/*
 * This program computes the distribution of a variable in a file 
 * it reads the indicated input file, extracting the indicated variable,
 * and sorts occurances into bins, and prints out the value in each
bin.
 *
 * Linas Vepstas November 1992 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef ANSI_C
#include <stdlib.h>
#else
#define EXIT_FAILURE 0
#endif

#define DEBUG 0

/* =================================================== */

#define ARRAY_SIZE 1000000

int n=0;
double xmin, xmax;
double x[ARRAY_SIZE];

/* =================================================== */
/* 
 * this subroutine scans the input file 
 */
#ifdef ANSI_C
void read_file (char *file_name, char *var_id)
#else
void read_file (file_name, var_id)
char *file_name;
char *var_id;
#endif
{
   FILE *infile;
   char data[512];
   char delimiters[100];
   char * token;

   /* open the input file */
   infile = fopen (file_name,"r");
   if (infile == NULL) {
      fprintf (stderr, "Error: Cannot open file %s \n", file_name);
      exit (EXIT_FAILURE);
   }

   /* define the token delimiters */
   strcpy(delimiters,"\n\t ");

   /* initialize array variables */
   n = 0;
   xmax = -1.0e30;
   xmin = 1.0e30;

   /* read until end-of-file */
   while (!feof (infile)) {

      /* get new line */
      fgets (data, 512, infile);

      /* parse the string */
      token = strtok (data, delimiters);

      /* comment card  --- skip */
      if (!strcmp (token,"#")) continue; 

      /* search for the desired token */
      while (strcmp (token, var_id)) {
         token = strtok (NULL, delimiters);
         if (!strcmp (token,"#")) { token = NULL;  break; }
         if (token == NULL) break;
      }

      /* didin't find the desired token */
      if (token == NULL) continue;

      /* ahhh - we did find the desired token */
      token = strtok (NULL, delimiters);
      if (token == NULL) {
         fprintf (stderr, "Error: file format error n = %d \n", n);
         continue;
      }

      x[n] = atof (token);
      if (x[n] > xmax) xmax = x[n];
      if (x[n] < xmin) xmin = x[n];

      n++;
      if (n > ARRAY_SIZE) {
         fprintf (stderr, "Error: array size exceeded \n");
         break;
      }
   }

   fprintf (stderr, "Done reading the input file \n");
   fclose (infile);
}


/* =================================================== */
/*
 * This routine sorts the data into bins
 */

#ifdef ANSI_C
void bin_sort (int nbins)
#else
void bin_sort (nbins)
int nbins;
#endif
{
   double bin_size;
   int i, bin;
   int * bins;
   double xbinmin, xbinmax;

   /* OK, in this section, we split the data up into bins */

   /* first, determine the number of bins required */
   xmax += 1.0e-8;
   xmin -= 1.0e-8;
   bin_size = (xmax - xmin) / ((double) nbins);

   printf ("# nbins = %d \n", nbins);
   printf ("# bin size = %f \n", bin_size);
   printf ("# num items sorted = %d \n", n);

   /* alloc the bins */
   bins = (int *) malloc (nbins * sizeof (int));

   /* initialize */
   for (i=0; i<nbins; i++) bins[i] = 0;

   /* sort */
   for (i=0; i<n; i++) {
      bin = (int) ((x[i] - xmin) / bin_size);
      bins[bin] ++;
   }

   /* print */
   for (i=0; i<nbins; i++) {
      xbinmin = xmin + ((double) i) * bin_size;
      xbinmax = xmin + ((double) i+1) * bin_size;
      printf (" %.3f %d \n", xbinmin, bins[i]);
   }

   fflush (stdout);

}

/* =================================================== */

#ifdef ANSI_C
main (int argc, char *argv[])
#else
main (argc, argv)
int argc;
char *argv[];
#endif
{
   if (argc != 4) {
      fprintf (stderr, "Usage: %s <infile> <var> <num bins> \n", 
            argv[0]);
      exit (EXIT_FAILURE);
   }


   printf ("# processing file %s \n", argv[1]);
   read_file (argv[1], argv[2]);

   printf ("# bin sorting for %s \n", argv[2]);
   bin_sort (atoi(argv[3]));

   exit (0);
}

/* ====================== END OF FILE ===================== */
