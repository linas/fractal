
#include "fileio.h"

#define XGRAPH

/* ========================================================== */
/* 
 * This routine acts as a filter for my data format, pulling out the
 * indicated data values, and dumping them into an xmgr nxy data format.
 *
 * Linas Vepstas November 1992
 */

#ifdef ANSI_C
void read_n_dump (char *file_name, char *xlabel, char *ylabel)
#else
void read_n_dump (file_name, xlabel, ylabel)
char *file_name, *xlabel, *ylabel;
#endif
{
   FILE *fi;
   char current_record[256];
   char delimiters[20];
   char *token;
   double x, y;

#ifdef XGRAPH
   /* print a data set label for xgraph */
   printf ("\"%s\n", ylabel);
#endif

   /* open the file */
   fi = fopen (file_name,"r");
   if (fi == NULL) {
      fprintf (stderr, "Error: uanble to open file %s \n", file_name);
      exit (EXIT_FAILURE);
   }

   /* set up delimiters for string parsing */
   strcpy(delimiters,"\n\t ");

   /* loop until end-of-file */
   while (fgets (current_record, 256, fi)) {

      /* comment card */
      if (current_record[0] == '#') continue;


      /* parse the record into tokens */
      token = strtok (current_record, delimiters);

      /* while not end-of-line */
      while (token) {

         /* comment card  --- skip */
         if (!strcmp (token,"#")) continue;

         /* does this token match the desired x label ?? */
         if (!strcmp (xlabel, token)) {

            /* if so, update x */
            token = strtok (NULL, delimiters);
            x = atof (token);
            continue;
         }

         /* does this token match the desired y label ?? */
         if (!strcmp (ylabel, token)) {

            /* if so, print out x and y */
            token = strtok (NULL, delimiters);
            y = atof (token);
            printf ("%f %f\n", x, y);
            continue;
         }

         token = strtok (NULL, delimiters);
      }
   }

   fclose (fi);

#ifdef XGRAPH
   /* print a blank line for xgraph */
   printf ("\n");
#endif
   fflush (stdout);

}


/* ========================================================== */
/* 
 * This routine acts opens and reads a file in my format, extracting the
 * indicated data.
 *
 * Linas Vepstas December 1992
 */

#define MSIZE 1000

#ifdef ANSI_C
int read_data (char *file_name, char *xlabel, char *ylabel, double **xdata, double **ydata)
#else
int read_data (file_name, xlabel, ylabel, xdata, ydata)
char *file_name, *xlabel, *ylabel;
double **xdata, **ydata;
#endif
{
   FILE *fi;
   char current_record[256];
   char delimiters[20];
   char *token;
   double x, y;
   double *xarr, *yarr;
   int nelts;
   int array_size;

   /* malloc initial data arrays */
   xarr = (double *) malloc (MSIZE * sizeof (double));
   yarr = (double *) malloc (MSIZE * sizeof (double));
   if ((xarr==NULL) || (yarr==NULL)) {
      fprintf (stderr, "read_data: Error: unable to malloc\n");
      exit (EXIT_FAILURE);
   }
   array_size = MSIZE;
   nelts = 0;

   /* open the file */
   fi = fopen (file_name,"r");
   if (fi == NULL) {
      fprintf (stderr, "Error: uanble to open file %s \n", file_name);
      exit (EXIT_FAILURE);
   }

   /* set up delimiters for string parsing */
   strcpy(delimiters,"\n\t ");

   /* loop until end-of-file */
   while (fgets (current_record, 256, fi)) {

      /* comment card */
      if (current_record[0] == '#') continue;

      /* parse the record into tokens */
      token = strtok (current_record, delimiters);

      /* while not end-of-line */
      while (token) {

         /* comment card  --- skip */
         if (!strcmp (token,"#")) continue;

         /* does this token match the desired x label ?? */
         if (!strcmp (xlabel, token)) {

            /* if so, update x */
            token = strtok (NULL, delimiters);
            x = atof (token);
            continue;
         }

         /* does this token match the desired y label ?? */
         if (!strcmp (ylabel, token)) {

            /* if so, store x and y */
            token = strtok (NULL, delimiters);
            y = atof (token);
            xarr[nelts] = x;
            yarr[nelts] = y;
            nelts ++;

            /* make sure we've got enough memory */
            if (nelts >= array_size) {
               array_size += MSIZE;
               xarr = (double *) realloc (xarr, array_size * sizeof (double));
               yarr = (double *) realloc (yarr, array_size * sizeof (double));

               if ((xarr==NULL) || (yarr==NULL)) {
                  fprintf (stderr, "read_data: Error: unable to malloc\n");
                  exit (EXIT_FAILURE);
               }
            }
            continue;
         }

         token = strtok (NULL, delimiters);
      }
   }

   fclose (fi);

   /* Be nice, and realloc to the final array size */
   xarr = (double *) realloc (xarr, nelts * sizeof (double));
   yarr = (double *) realloc (yarr, nelts * sizeof (double));

   if ((xarr==NULL) || (yarr==NULL)) {
      fprintf (stderr, "read_data: Error: unable to malloc\n");
      exit (EXIT_FAILURE);
   }

   *xdata = xarr;
   *ydata = yarr;
   
   return (nelts);
}

/* ========================================================== */
/* 
 * This routine searches a file for the first occurance of "token"
 * and returns the associated value.
 *
 * Linas Vepstas November 1992
 * modified December 1993
 */

#ifdef ANSI_C
int read_value (char *file_name, char *value_name, double *value)
#else
int read_value (file_name, value_name, value)
char *file_name, *value_name;
double *value;
#endif
{
   FILE *fi;
   char current_record[256];
   char delimiters[20];
   char *token;

   /* open the file */
   fi = fopen (file_name,"r");
   if (fi == NULL) {
      fprintf (stderr, "Error: uanble to open file %s \n", file_name);
      exit (EXIT_FAILURE);
   }

   /* set up delimiters for string parsing */
   strcpy(delimiters,"\n\t ");

   /* loop until end-of-file */
   while (fgets (current_record, 256, fi)) {

      /* comment card */
      if (current_record[0] == '#') continue;

      /* parse the record into tokens */
      token = strtok (current_record, delimiters);

      /* while not end-of-line */
      while (token) {

         /* comment card  --- skip */
         if (!strcmp (token,"#")) continue;

         /* does this token match the desired label ?? */
         if (!strcmp (value_name, token)) {

            /* if so, update x */
            token = strtok (NULL, delimiters);
            *value = atof (token);
            fclose (fi);
            /* value found, return sucess */
            return (0);
         }

         token = strtok (NULL, delimiters);
      }
   }

   /* value not found, return error */
   return (-1);
}

/* ====================== END OF FILE ======================= */
