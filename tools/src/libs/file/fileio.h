

#ifdef AIX221
#ifndef _C_func
#define _C_func
#endif
#define EXIT_SUCESS 0
#define EXIT_FAILURE 1
#endif

#ifdef AIX315
#define ANSI_C
#endif

#ifdef ANSI_C
#include <stdlib.h>
#endif


#include <string.h>
#include <stdio.h>
#include <math.h>

/* ========================================================== */
/* function prototypes */

#ifdef ANSI_C

extern void read_n_dump (char *filename, char *xlabel, char *ylabel);
extern int read_data (char *filename, char *xlabel, char *ylabel, double **xdata, double **ydata);
extern int read_value (char * file_name, char *value_name, double *value);

#else

extern void read_n_dump ();
extern int read_data ();
extern int read_value ();

#endif


/* ====================== END OF FILE ======================= */
