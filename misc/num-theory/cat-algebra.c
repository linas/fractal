
/*
 * cat-algebra.c
 *
 * Categorical algebra mobius to num-theo mobius
 *
 * Linas Vepstas November 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

	double fk = moebius_mu(k);

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stdout, "Usage: %s <nbins> <kmax>\n", argv[0]);
		exit(1);
	}
	int nbins = atoi (argv[1]);
	int kmax = atoi(argv[2]);

	rebin (nbins, kmax);
}
