
#include "fileio.h"

#define XGRAPH

/* ========================================================== */
/* 
 * This program acts as a filter for my data format, pulling out the
 * indicated data values, and dumping them into an xmgr nxy data format.
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
   int i;

   /* Usage */
   if (argc < 4) {
      fprintf (stderr, 
          "Usage: %s <filename> <x-axis label> <y-axis label> [<y axis> ...] \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   /* get the filename */
   file_name = argv[1];
   xlabel = argv[2];

   for (i=3; i<argc; i++) {
      ylabel = argv[i];
      read_n_dump (file_name, xlabel, ylabel);
   }

   exit (0);
}
/* ====================== END OF FILE ===================== */
