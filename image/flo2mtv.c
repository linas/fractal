/* 
 * FUNCTION:
 * Convert .flo files to .mtv files
 * Also convert .pfm to .ppm files.
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 * Fix colormap 16 Dec 2017
 * Add pfm support January 2024
 */

#include <endian.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct rgb {
   char r;
   char g;
   char b;
};

static struct rgb vlt[256];

/* ------------------------------------------------------------ */

void make_cmap (void)
{
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
            // buggy, but used for everything before 2018
            // vlt[i].g = (char) (510 - 2*i);
            vlt[i].g = (char) (525 - 2*i);
            vlt[i].b = (char) (i-180)/3;
        }
}

/*-------------------------------------------------------------------*/
// My version of fgets, because fgets is sometimes broken on other machines.
// I don't remember what the problem was.
static char *my_fgets(char *s, int size, FILE *stream)
{
   for (int i=0; i< size; i++)
   {
      int ch =  fgetc(stream);
      s[i] = ch;
      if ('\n' == ch) break;
      if (0 == ch) break;
   }
   return s;
}

int main (int argc, char* argv[])
{
   if (argc > 2) {
      printf ("Usage: %s <flo filename> \n", argv[0]);
      exit (1);
   }

   /* open up the file */
   FILE *filh;
   if (argc == 2) {
      filh = fopen (argv[1], "r");
   } else {
      filh = stdin;
   }

   /* read file magic */
   char ascii_string[80];
	my_fgets(ascii_string, 80, filh);

#define FLO_FILE 1
#define PFM_FILE 2
	int ftype = FLO_FILE;
	if ('P' == ascii_string[0])
		ftype = PFM_FILE;

   int width, height;
	if (FLO_FILE == ftype)
	{
		sscanf (ascii_string, "%d %d", &width, &height);
	}
	else if (PFM_FILE == ftype)
	{
		my_fgets(ascii_string, 80, filh);
		sscanf (ascii_string, "%d %d", &width, &height);

		// Ignore the "scale factor"
		my_fgets(ascii_string, 80, filh);
	}
	else
	{
		fprintf(stderr, "%s: unknown file type\n", argv[0]);
		exit(1);
	}

   make_cmap ();

	if (FLO_FILE == ftype)
	{
		// Write an mtv file, which is just width, height and data
		fprintf (stdout, "%d %d\n", width, height);
	}
	else if (PFM_FILE == ftype)
	{
		// Wrte a ppm file
		fprintf (stdout, "P6\n%d %d\n255\n", width, height);
	}
   
   float *in_row = (float *) malloc (width*sizeof (float));
   char *out_row = (char *) malloc (3*width*sizeof (char));

   /* loop over pixels, using bogus colormap */
   for (int i=0; i<height; i++) {
      fread (in_row, sizeof(float), width, filh);
      for (int j=0; j<width; j++) {
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
