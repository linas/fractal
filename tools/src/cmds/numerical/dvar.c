
#include "fileio.h"
#include "scale.h"
#include <math.h>

#define XGRAPH
#define NSPEC 300

/* ========================================================== */
/* 
 * This program reads the indicated file, computes the spectral density
 * of the data, and dumps the spectrum to stdout.
 *
 * Linas Vepstas November 1992
 */

#ifdef ANSI_C
int main (int argc, char *argv[])
#else
int main (argc, argv)
int argc;
char *argv[];
#endif
{
   char *file_name, *xlabel, *ylabel;
   int i, j;
   double *xdata, *ydata;
   int npoints;
   double sigma, S;
   double sigmaZero, sigmaMax;
   double logsigma;
   double scale;

   /* Usage */
   if (argc < 4) {
      fprintf (stderr, 
          "Usage: %s <filename> <x-axis label> <y-axis label> [<y axis> ...] \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   /* get the filename */
   file_name = argv[1];
   xlabel = argv[2];

   /* loop over all files */
   for (i=3; i<argc; i++) {
      ylabel = argv[i];

#ifdef XGRAPH
      /* print a data set label for xgraph */
      printf ("\"log10 %s\n", ylabel);
#endif

      npoints = read_data (file_name, xlabel, ylabel, &xdata, &ydata);

      /* The smallest interesting sigma is that which characterizes the
       * difference between neighboring data points -- i.e. sigma = 1.0. 
       * The largest interesting sigma is that which spans most of the
       * data set -- say, ahh, one-twentieth of the total number of data 
       * points. (after that, you just get bad data) */
      sigmaZero = 1.0;
      sigmaMax = 0.05 * ((double) npoints);

      /* to get evenly spaced points on a log graph, thake the ratio of
       * the endpoints of the graph, and take nth root, where n is the 
       * number of graph points. */
      scale = pow ((sigmaMax / sigmaZero), (1.0/((double) NSPEC)));

      sigma = sigmaZero;

      /* compute derivative varience  */
      for (j=0; j<NSPEC; j++) {
         sigma *= scale;
         S = DifferentialRMS (ydata, npoints, sigma); 
         S = log10 (S);
         logsigma = log10 (sigma);

         printf ("%f %f \n", logsigma, S);
         fflush (stdout);
      }

#ifdef XGRAPH
      /* print a blank line for xgraph */
      printf ("\n");
#endif
      fflush (stdout);

   }
   exit (0);
}
/* ====================== END OF FILE ===================== */
