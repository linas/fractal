
/* do stuff */

#include <math.h>
#include <stdio.h>

#include "bag_ener.h"

int main ()
{
	int i;
	double theta;
	int ispect;
	int k;
	int kpty;
	double emax;
	int nfound;

	#define NLVLS 100
	double energy_levels[NLVLS];

	emax = NLVLS * 3;

	theta = 0.1;
	ispect = +1;
	k = 0;
	nfound = quark_energy (energy_levels, theta, ispect, k, emax, kpty);

	printf ("found %d levels: \n", nfound);
	for (i=0; i<nfound; i++)
	{
		printf ("%d %g\n", i, energy_levels[i]);
	}
}
