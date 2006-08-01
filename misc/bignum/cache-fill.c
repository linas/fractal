/*
 * cache-fill.c
 *
 * File the file database with high-precision zeta values
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mp_zeta.h"

/* ==================================================================== */

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [start]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* place to start */
	int nstart = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* And pad a bit, just to be safe... */
	int bits = (int) (v + 1100);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_t term;
	mpf_init (term);

	int n;
	// for (n=nstart; n>1 ; n--)
	for (n=2; n<nstart; n++)
	{
		time_t start = time(0);
		fp_zeta (term, n, prec);
		time_t end = time(0);
		int elapsed = end-start;
		
		printf ("%d\t", n);
		mpf_sub_ui (term, term, 1);
		fp_prt ("", term);
		printf ("\t%d\n",elapsed);
		fflush (stdout);
	}
	
	return 0;
}

