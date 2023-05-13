
/* do stuff */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bessel.h"
#include "bag_ener.h"

int main (int argc, char * argv[])
{
	int i;
	double theta;
	int k;
	double emax;
	int nfound;

	#define NLVLS 1000
	double pos_even_levels[NLVLS];
	double pos_odd_levels[NLVLS];
	double neg_even_levels[NLVLS];
	double neg_odd_levels[NLVLS];

	emax = NLVLS * 3;

	if (3 > argc)
	{
		printf ("Usage: %s theta emax\n", argv[0]);
		exit (1);
	}
	theta = atof (argv[1]);
	emax = atof (argv[2]);


	printf ("#\n");
	printf ("# FILE: \n");
	printf ("#\n");
	printf ("# even levels K=L      kpty = +1 \n");
	printf ("# odd  levels K=L+/-1  kpty = -1 \n");
	printf ("theta %21.15g\n", theta);
	printf ("emax %21.15g\n", emax);

	for (k=0; k<emax; k++)
	{
		int nfoundp, nfoundn, np, nn;

		np = quark_energy (pos_even_levels, theta, 1, k, emax, 1);
		nn = quark_energy (pos_odd_levels, theta, 1, k, emax, -1);
		nfoundp = fmin (np, nn);
		np = quark_energy (neg_even_levels, theta, -1, k, emax, 1);
		nn = quark_energy (neg_odd_levels, theta, -1, k, emax, -1);
		nfoundn = fmin (np, nn);
		nfound = fmin (nfoundp, nfoundn);

		printf ("k %d nlevels %d\n", k, nfound);
		for (i=0; i<nfound; i++)
		{
			printf ("%d\t%21.15g\t%21.15g\t%21.15g\t%21.15g\n", 
				i, 
				pos_even_levels[i], pos_odd_levels[i],
				neg_even_levels[i], neg_odd_levels[i]);
		}
	}
	return 0;
}
