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
#include "unutil.c"

double fnorm(double beta)
{
	double sum = 0.0;
	double bei = 1.0;
	for (int i=0; i<1000; i++)
	{
		sum += t_k(beta, i) / bei;
		bei *= beta;
		if (1.0 < bei * 1e-17) break;
	}
	return sum;
}

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

double fint(double beta)
{
#define FSAMP 1513
	double sum = 0.0;
	for (int i=0; i< FSAMP; i++)
	{
		double x = (((double) i) + 0.5) / ((double) FSAMP);
		double boo = 1.0+x;
		if (beta < boo) break;

		double fno = fnorm(boo);
		if (0.0 < fno)
			sum += 1.0/fno;
	}
	return sum / ((double) FSAMP);
}

int main(int argc, char* argv[])
{
	int imax = 800;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
		double eff = 1.0 / fnorm(beta);
		double inf = fint(beta);
		printf("%d	%f	%f	%f", i, beta, eff, inf);
		// for (int n=0; n<6; n++)
		for (int n=0; n<19; n+=3)
		{
			double eno = (n+1) * enorm(beta, n);
			eno = 1.0 / eno;
			printf("	%f", eno);
		}

		printf("\n");
		fflush(stdout);
	}
}
