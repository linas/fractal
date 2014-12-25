
/**
 * fourier.c
 *
 * Harmonic analysis of the a_k series
 */

#include <math.h>
#include <stdio.h>

#include "falling.h"

double four(double omega)
{
	int k;

	double sum = 0.0;
	for (k=1; k< 100; k++)
	{
		double ak = a_k(k);
		double r = exp(-3.3*sqrt(k));
		ak *= r;

		sum += ak * sin(omega*k);
	}
	return sum;
}

int main(int argc, char* argv[])
{
	double om;

	for (om=2.8; om < 3.1; om += 0.01)
	{
		printf(" its %g %g\n", om, four(om));
	}
}
