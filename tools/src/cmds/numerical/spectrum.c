
#include "fileio.h"
#include "spectra.h"
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
   double omega, S;
   double omegaMax, omegaZero;
   double scale;
   double logfreq;

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
      printf ("\"%s\n", ylabel);
#endif

      npoints = read_data (file_name, xlabel, ylabel, &xdata, &ydata);

      /* The lowest interesting frequency is one-quarter wave spread
       * over all of the data points -- i.e. pi-halves divided by number
       * of data points.  The highest interesting frequency is the 
       * (Nyquist??) frequency, which is data points on half-cycles 
       * -- i.e. pi.  Above this frequency, aliasing maps higher 
       * frequencies into the lower band. */
      omegaZero = M_PI / (2.0 * ((double) npoints));
      omegaMax = M_PI;

      /* to get evenly spaced points on a log graph, thake the ratio of
       * the endpoints of the graph, and take nth root, where n is the 
       * number of graph points. */
      scale = pow ((omegaMax / omegaZero), (1.0/((double) NSPEC)));

      omega = omegaZero;
      /* compute spectral density of each */
      for (j=1; j<NSPEC; j++) {
         omega *= scale;
         /* S = SpectralDensity (ydata, npoints, omega); */
         S = SmoothSpectralDensity (ydata, npoints, omega);

         S = log10 (S);
         logfreq = log10 (omega);

         printf ("%f %f \n", logfreq, S);
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
