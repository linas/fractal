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

double en0(double beta)
{
#define NSAMP 3512
	double sum = 0.0;
	for (int i=0; i< NSAMP; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NSAMP);
		sum += 1.0 / gp_invar(beta, x);
	}
	return sum / ((double) NSAMP);
}

int main(int argc, char* argv[])
{
	int imax = 514;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
		double no = en0(beta);
		printf("%d	%f	%f\n", i, beta, no);
		fflush(stdout);
	}
}
