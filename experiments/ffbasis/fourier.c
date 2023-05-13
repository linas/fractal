
/**
 * fourier.c
 *
 * Harmonic analysis of the a_k series
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "falling.h"

double four(double omega, double phi)
{
	int k;

	double sum = 0.0;
	for (k=200; k< 300; k++)
	{
		double ak = a_k(k);
		double r = exp(-3.3*sqrt(k));
		ak *= r;

		sum += ak * sin(omega*k + phi);
	}
	return sum;
}

int main(int argc, char* argv[])
{
	double om, phi;

	// phi = atof(argv[1]);

	// for (om=2.88; om < 2.91; om += 0.0002)
	// for (om=2.85; om < 2.95; om += 0.0002)
	for (om=3.0; om < 3.13; om += 0.0002)
	{
		printf("%g", om);
		for (phi=0; phi<6.283; phi += 1.0) {
			printf("	%g", four(om, phi));
		}
		printf("\n");
	}
}
