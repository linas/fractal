
#include "fileio.h"
#include "var.h"
#include <math.h>

#define XGRAPH

#define EXIT_SUCCESS 0
#define MIN(x,y) ((x<y)?y:x);

/* ========================================================== */
/* 
 * This program reads the indicated file, and computes the correlation
 * coefficient between the two indicated values.  The correlation is
 * computed with a sliding gaussian filter, so that we can check whther the
 * correlation for the time series is about constant, or whether it changes
 * dramitically.
 *
 * Linas Vepstas Decemebr 1993
 */

#ifdef ANSI_C
int main (int argc, char *argv[])
#else
int main (argc, argv)
int argc;
char *argv[];
#endif
{
   char *file_name, *xlabel;
   char *corr_a_label, *corr_b_label;
   double *xdata, *adata, *bdata;
   int npoints_a, npoints_b, npoints;
   int i, j;
   double r;
   double width;

   /* Usage */
   if ((argc < 6) || (argc %2 != 0)) {
      fprintf (stderr, 
          "Usage: %s <filename> <width> <x-axis label> <corr-a label> <corr-b label> [...] \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   /* get the filename */
   file_name = argv[1];
   width = atof (argv[2]);
   xlabel = argv[3];

   /* loop over all files */
   for (i=4; i<argc; i+=2) {
      corr_a_label = argv[i];
      corr_b_label = argv[i+1];

#ifdef XGRAPH
      /* print a data set label for xgraph */
      printf ("\"corr <%s %s>\n", corr_a_label, corr_b_label);
#endif

      npoints_a = read_data (file_name, xlabel, corr_a_label, &xdata, &adata);
      free (xdata);
      npoints_b = read_data (file_name, xlabel, corr_b_label, &xdata, &bdata);
      npoints = MIN (npoints_a, npoints_b);


      /* loop over data */
      for (j=0; j<npoints; j++) {
         r = GCorrelation (adata, bdata, npoints, (double) j, width);
         printf ("%f %f \n", xdata [j], r);
         fflush (stdout);
      }

#ifdef XGRAPH
      /* print a blank line for xgraph */
      printf ("\n");
#endif
      fflush (stdout);
      free (xdata);
      free (adata);
      free (bdata);

   }
   exit (EXIT_SUCCESS);

}

/* ====================== END OF FILE ===================== */
