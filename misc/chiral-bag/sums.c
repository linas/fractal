
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

	double t;
	double barn, k0_barn;
	double theo;

	#define NLVLS 1000
	double pos_even_levels[NLVLS];
	double pos_odd_levels[NLVLS];
	double neg_even_levels[NLVLS];
	double neg_odd_levels[NLVLS];

	emax = NLVLS * 3;

	theta = 0.3;
   emax = 600;

	if (2 > argc)
	{
		printf ("Usage: %s theta\n", argv[0]);
		exit (1);
	}
	theta = atof (argv[1]);

	theo = (theta - sin(theta)*cos(theta)) / M_PI;
	printf ("theta=%g theo=%g\n", theta, theo);

	for (t=0.5; t; t *= sqrt(0.5))
	{
		emax = sqrt (- log (1.0e-17)) / t;
		if (emax > 3*NLVLS) break;
	
		barn = 0.0;
		k0_barn = 0.0;
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
	
			for (i=0; i<nfound; i++)
			{
				double en;
				en = pos_even_levels[i];
				barn += exp (- t*t*en*en);
				en = pos_odd_levels[i];
				barn += exp (- t*t*en*en);
				en = neg_even_levels[i];
				barn -= exp (- t*t*en*en);
				en = neg_odd_levels[i];
				barn -= exp (- t*t*en*en);
			}
			if (k==0) k0_barn = barn;
		}
		barn *= -0.5;
		k0_barn *= -0.5;

		printf ("t=%8.5g	emax=%g	bar=%g	k0=%g diff=%g\n", 
			t, emax, barn, k0_barn, barn-k0_barn);
	}
	return 0;
}
