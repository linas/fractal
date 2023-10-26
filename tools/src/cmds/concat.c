

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


   if (argc  < 3) {
      fprintf (stderr, "%s concatenates file 1 and file 2 \n", argv[0]);
      fprintf (stderr, "Usage: %s <filename 1> <filename 2> \n", argv[0]);
      exit (EXIT_FAILURE);
   }

   for (i=1; i<argc; i++) {

      file_name = argv [i];
      
      /* open the file */
      fi = fopen (file_name,"r");
      if (fi == NULL) {
         fprintf (stderr, "Error: uanble to open file %s \n", file_name);
         fflush (stdout);
         exit (EXIT_FAILURE);
      }

      while (!feof (fi)) {
         c = fgetc (fi);
         putc (c, stdout);
      }

      fclose (fi);

   }
   fflush (stdout);

}

/* ====================== END OF FILE ======================= */
