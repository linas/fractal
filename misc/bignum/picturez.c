/* 
 * picturez.c
 *
 * Crazy pictures of digits of zeta.
 *
 * Linas Vepstas July 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#include "mp_zeta.h"

int
main (int argc, char * argv[])
{
	if (argc<4)
	{
		fprintf (stderr,"Usage: %s <outfile> <pixel width> <pixel height>\n", argv[0]);
		exit (1);
	}
	char * outfile = argv[1];
	int width = atoi (argv[2]);
	int height = atoi (argv[3]);
	int base = 10;

	/* Compute number of binary bits this corresponds to. */
	int prec = width;
	double v = ((double) prec) *log(10.0) / log(2.0);
	int bits = (int) (v + 100);
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);

	float *arr = (float *) malloc (width*height*sizeof(float));
	
	mpf_t val;
	mpf_init (val);
	int n;
	for (n=2; n<= height+1; n++)
	{
		char buff[2000];
		mp_exp_t ep = 0;
		fp_zeta (val, n, width+50);
		mpf_get_str (buff, &ep, base, width+10, val);
		int i=n-2;
		int j;
		for (j=0;j<width; j++)
		{
			arr[i*width+j] = ((float) (buff[j]-'0')) / ((float) base);
		}
	}

	/* dump the floating point data */
	FILE * fp = fopen (outfile, "w");
	fprintf (fp, "%d %d\n", width, height);
	fwrite (arr, sizeof(float), width*height, fp);
	fclose (fp);

	return 0;
}
