/*
 * ortho.c
 *
 * Orthonormality experiments.
 * January 2024
 */

#include <math.h> 
#include <stdio.h>
#include <stdlib.h> 

#include "unref.c"

// Normalization
double enorm(double beta, int n)
{
#define NSAMP 6513
	double sum = 0.0;
	for (int i=0; i< NSAMP; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NSAMP);
		double xn = pow(x, n);
		sum += xn / gp_invar(beta, x);
	}
	return sum / ((double) NSAMP);
}

int main(int argc, char* argv[])
{
	int imax = 800;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
		printf("%d	%f", i, beta);
		for (int n=0; n<6; n++)
			printf("	%f", enorm(beta, n));

		printf("\n");
		fflush(stdout);
	}
}
