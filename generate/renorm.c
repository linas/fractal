/*
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "opers.h"

extern FILE *Fopen();
extern FILE *Fopenr();

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

/*-------------------------------------------------------------------*/
// Read a *.flo greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
float* read_flo_file(const char *fname,
                     unsigned int* data_width, unsigned int* data_height)
{
	/* open input file */
	FILE *fh = Fopenr (fname, ".flo");
	if (NULL == fh) return NULL;

#define CLEN 80
   char str[CLEN];
	my_fgets(str, CLEN, fh);

	sscanf(str, "%d %d", data_width, data_height);
	fprintf(stderr, "Input file %s has width %d height %d \n",
	        fname, *data_width, *data_height);

	size_t datalen = (*data_width) * (*data_height);
	float* data_in = (float *) malloc (datalen * sizeof(float));

	/* read floating point data */
	fread(data_in, sizeof(float), datalen, fh);
	fclose(fh);

	return data_in;
}

/*-------------------------------------------------------------------*/
// Read a *.pfm greyscale floating-point pixmap.
// This reads the PFM file format as described in various places.
float* read_pfm_file(const char *fname,
                     unsigned int* data_width, unsigned int* data_height)
{
	/* open input file */
	FILE *fh = Fopenr (fname, ".pfm");
	if (NULL == fh) return NULL;

#define PFLEN 80
   char str[PFLEN];
	my_fgets(str, PFLEN, fh);

	// Handle only greyscale
	if ('P' != str[0] || 'f' != str[1]) return NULL;

	// Read the dimensions
	my_fgets(str, PFLEN, fh);
	sscanf(str, "%d %d", data_width, data_height);
	fprintf(stderr, "Input file %s has width %d height %d \n",
	        fname, *data_width, *data_height);

	// Read the "scale factor" and ignore it.
	my_fgets(str, PFLEN, fh);
	// double scale = atof(str);

	size_t datalen = (*data_width) * (*data_height);
	float* data_in = (float *) malloc (datalen * sizeof(float));

	/* read floating point data */
	fread(data_in, sizeof(float), datalen, fh);
	fclose(fh);

	return data_in;
}

/*-------------------------------------------------------------------*/
// Write a *.flo greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
void write_flo_file(const char *fname,
                    float* data,
                    unsigned int data_width, unsigned int data_height)
{
   FILE *fh = Fopen(fname, ".flo");
   if (NULL == fh)
	{
      fprintf (stderr, "Can't open output file %s.flo\n", fname);
		return;
   }

   /* dump the size to output */
   fprintf (fh, "%d %d\n", data_width, data_height);
   fwrite (data, sizeof(float), data_width*data_height, fh);
   fclose (fh);
}

/*-------------------------------------------------------------------*/
// Write a *.pfm greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
void write_pfm_file(const char *fname,
                    float* data,
                    unsigned int data_width, unsigned int data_height)
{
	FILE *fh = Fopen(fname, ".pfm");
	if (NULL == fh)
	{
		fprintf (stderr, "Can't open output file %s.pfm\n", fname);
		return;
	}

	/* dump the size to output */
	fprintf (fh, "Pf\n");
	fprintf (fh, "%d %d\n", data_width, data_height);
	double scale = 1.0;
	// if (htonl(47) != 47) scale = -1.0; // Little-endian
	if (htobe16(47) != 47) scale = -1.0; // Little-endian
	fprintf (fh, "%f\n", scale);

	/* And now the data. */
	fwrite (data, sizeof(float), data_width*data_height, fh);
	fclose (fh);
}

/*-------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
   unsigned int data_width, data_height;/* data array dimensions */
   float 		scale_fact;
   float 		offset;

   if (argc < 3) {
      fprintf (stderr, "Usage: %s <input file> <output file> [<scale factor>] [<offset>] \n", argv[0]);
      return 1;
   }

   fprintf (stderr, "%s %s %s ", argv[0], argv[1], argv[2]);

   if (4 <= argc) {
      scale_fact = (float) atof (argv[3]);
      fprintf (stderr, "%g ", scale_fact);
   } else {
      scale_fact = 1.0;
   }

   if (5 <= argc) {
      offset = (float) atof (argv[4]);
      fprintf (stderr, "%g ", offset);
   } else {
      offset = 0.0;
   }
   fprintf (stderr, "\n");

	/*-----------------------------------------------*/
	/* open input file */
	const char* suff = ".flo";
	float* data_in = read_flo_file(argv[1], &data_width, &data_height);
	if (NULL == data_in)
	{
		suff = ".pfm";
		data_in = read_pfm_file(argv[1], &data_width, &data_height);
	}
	if (NULL == data_in)
	{
		fprintf (stderr, " Can't open input file %s \n", argv[1]);
		return 2;
	}

   /*-----------------------------------------------*/
	/* Strip out the path */
	char * p = strdup(argv[0]);
	char * name = strrchr(p, '/');
	if (NULL == name) name = p;
	else name++;

   if (!strcmp ("takelog", name)) takelog (data_in, data_width, data_height);
   else if (!strcmp ("recip", name)) recip (data_in, data_width, data_height);
   else if (!strcmp ("abs", name)) absolval (data_in, data_width, data_height);
   else  if (!strcmp ("takeroot", name)) takeroot (data_in, data_width, data_height, scale_fact);
   else  if (!strcmp ("renorm", name)) expand (data_in, data_width, data_height, scale_fact, offset);
   else  if (!strcmp ("reclamp", name)) rescale (data_in, data_width, data_height, scale_fact);
   else  if (!strcmp ("clamp", name)) clamp (data_in, data_width, data_height, scale_fact, offset);
   else  if (!strcmp ("dump", name)) dump (data_in, data_width, data_height);
	else
	{
		fprintf(stderr, "Error: unknown transform %s\n", name);
		exit(1);
	}

/*
   aver = avg (data_in, data_width, data_height);
   msq = sqdev (data_in, data_width, data_height);
   msq = sqrt (msq);
   printf ("yo average was %f rms %f \n", aver, msq);
   expand (data_in, data_width, data_height, scale_fact/msq, aver);
*/

	if (0 == strcmp(suff, ".flo"))
		write_flo_file(argv[2], data_in, data_width, data_height);

	if (0 == strcmp(suff, ".pfm"))
		write_pfm_file(argv[2], data_in, data_width, data_height);

	return 0;
}

/* --------------- END OF FILE ------------------------- */
