
/* do stuff */

#include <math.h>
#include <stdio.h>

#include "bessel.h"
#include "bag_ener.h"

int main ()
{
	int i;
	double theta;
	int k;
	double emax;
	int nfound;

	double t;
	double baryon_number;
	double theo;

	#define NLVLS 1000
	double pos_even_levels[NLVLS];
	double pos_odd_levels[NLVLS];
	double neg_even_levels[NLVLS];
	double neg_odd_levels[NLVLS];

	emax = NLVLS * 3;

	theta = 1.0;
   emax = 600;

	for (t=0.5; t; t *= 0.5)
	{
		emax = sqrt (- log (1.0e-15)) / t;
		if (emax > 3*NLVLS) break;
	
		baryon_number = 0.0;
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
				baryon_number += exp (- t*t*en*en);
				en = pos_odd_levels[i];
				baryon_number += exp (- t*t*en*en);
				en = neg_even_levels[i];
				baryon_number -= exp (- t*t*en*en);
				en = neg_odd_levels[i];
				baryon_number -= exp (- t*t*en*en);
			}
		}
		baryon_number *= -0.5;
	
		theo = (theta - sin(theta)*cos(theta)) / M_PI;
		printf ("theta=%g    t=%g    bar=%g  theo=%g\n", 
			theta, t, baryon_number, theo);
	}
	return 0;
}
