
/* 
 * FUNCTION:
 * Convert .flo files to .mtv files
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <stdio.h>
#include <math.h>

struct rgb {
   char r;
   char g;
   char b;
};

static struct rgb vlt[256];

/* ------------------------------------------------------------ */

void make_cmap (void) {
    int i, j;
    struct rgb black;
    black.r = black.g = black.b = 0x0;
    for (i=0; i<256; i++) vlt[i] = black;

    /* set up a default look up table */
    /* ramp up to blue */
    for (i=0; i<60; i++) {
            vlt[i].r = 0;
            vlt[i].g = 0;
            vlt[i].b = (char) i*3;
        }
    /* ramp down from blue, up to green */
    for (i=60; i<120; i++) {
            vlt[i].r = 0;
            vlt[i].g = (char) (i-60)*3;
            vlt[i].b = (char) (120-i)*3;
        }
    /* ramp from green to yellow */
    for (i=120; i<180; i++) {
            /* vlt[i].r = (char) (((i-120)*7) / 2); */
            vlt[i].r = (char) (210 - (7*(180-i)*(180-i)) / 120);
            vlt[i].g = (char) (210 -i/4);
            vlt[i].b = 0;
        }
    /* ramp from yellow to red (pink) */
    for (i=180; i<240; i++) {
            vlt[i].r = (char) (210 + (3*(i-180))/4);
            vlt[i].g = (char) (510 - 2*i);
            vlt[i].b = (char) (i-180)/3;
        }
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   FILE *fil;
   int i, j, width, height;
   float *in_row;
   char *out_row;
   char ascii_string[80];

   if (argc > 2) {
      printf ("Usage: %s <flo filename> \n", argv[0]);
      exit (1);
   }

   make_cmap ();

   /* open up the file */
   if (argc == 2) {
      fil = fopen (argv[1], "r");
   } else {
      fil = stdin;
   }

   /* read the size */
   for (i=0; i<80; i++) {
      fread (&ascii_string[i], sizeof (char), 1, fil);
      if (ascii_string[i] == 0x0) break;
      if (ascii_string[i] == 0xa) break;
   }
   ascii_string[i] = 0x0; 

   sscanf (ascii_string, "%d %d", &width, &height);

   /* dump the size */
   fprintf (stdout, "%d %d\n", width, height);

   /* write null terminated string */
   /*  { char zip=0xa; fwrite (zip, sizeof(char), 1, stdout); }
   fflush (stdout); */
   
   in_row = (float *) malloc (width*sizeof (float));
   out_row = (char *) malloc (3*width*sizeof (char));

   /* loop over pixels, using bogus colormap */
   for (i=0; i<height; i++) {
      fread (in_row, sizeof(float), width, fil);
      for (j=0; j<width; j++) {
         int k;

         /* make sure large overflows don't wrap the integers */
         if (0.0 > in_row[j]) in_row[j] = 0.0;
         if (1.0 < in_row[j]) in_row[j] = 1.0;
         k = (int) (239.5 * in_row[j]);

         if (k<0) k=0;
         if (k>239) k=239;

         out_row[3*j] = vlt[k].r;
         out_row[3*j+1] = vlt[k].g;
         out_row[3*j+2] = vlt[k].b;
      }
      fwrite (out_row, 3*sizeof(char), width, stdout);
      fflush (stdout);
   }

   free (in_row);
   free (out_row);
}

/* ---------------------- END OF FILE ------------------------- */
