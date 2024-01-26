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

// Hausdorff moment (Hamburger moment on unit interval)
double hmoment(double beta, int n)
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
	int imax = 1200;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
		double eff = 1.0 / fnorm(beta);
		double inf = fint(beta);
		printf("%d	%f	%f	%f", i, beta, eff, inf);
		for (int n=0; n<6; n++)
		{
			double hmom = (n+1) * hmoment(beta, n);
			hmom = 1.0 / hmom;
			printf("	%f", hmom);
		}
		for (int n=6; n<19; n+=3)
		{
			double hmom = (n+1) * hmoment(beta, n);
			hmom = 1.0 / hmom;
			printf("	%f", hmom);
		}

		printf("\n");
		fflush(stdout);
	}
}
