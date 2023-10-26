
#define EXIT_SUCESS 0
#define EXIT_FAILURE 1

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ========================================================== */
/* 
 *
 * Linas Vepstas January 1993
 */

int main (int argc, char **argv)
{
   FILE *fi;
   char *file_name;
   char c;
   int i, N;


   if (argc != 3) {
      fprintf (stderr, "%s prints the first n characters of a file to standard out \n", argv[0]);
      fprintf (stderr, "Usage: %s <n> <filename> \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   N = atoi (argv[1]);
   file_name = argv[2];
      
   /* open the file */
   fi = fopen (file_name,"r");
   if (fi == NULL) {
      fprintf (stderr, "Error: uanble to open file %s \n", file_name);
      exit (EXIT_FAILURE);
   }

   /* if N is positive, print only the first N characters */
   /* if N is negative, print all but the first N characters */
   if (N>0) {
      for (i=0; i<N; i++) {
         c = fgetc (fi);
         putc (c, stdout);
      }
   } else {
      N = -N;
      for (i=0; i<N; i++) {
         c = fgetc (fi);
      }
      while (!feof (fi)) {
         c = fgetc (fi);
         putc (c, stdout);
      }
   }


   fclose (fi);

   fflush (stdout);

}

/* ====================== END OF FILE ======================= */
