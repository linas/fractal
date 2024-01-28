/*
 * moments.c
 *
 * Moments.
 * January 2024
 */

#include <math.h> 
#include <stdio.h>
#include <stdlib.h> 

// Hausdorff moment (Hamburger moment on unit interval)
double fmoment(double beta, int n)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow (tk, n);
		sum += tkn / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

int main(int argc, char* argv[])
{
	int imax = 1200;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
		printf("%d	%f", i, beta);
		for (int n=0; n<8; n++)
		{
			double fmom = fmoment(beta, n);
			printf("	%f", fmom);
		}
		printf("\n");
		fflush(stdout);
	}
}
