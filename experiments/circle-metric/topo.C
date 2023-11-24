/*
 * quick hack to create figures for the sigma algebra discussion.
 * November 2023
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	double omega = atof(argv[1]);
	double K = atof(argv[2]);

#define NPTS 2500
	for (int i=0; i<NPTS; i++)
	{
		double x = (i+0.5) / NPTS;
		double cx = x + omega - K * sin(2 * M_PI * x);
		cx -= floor(cx);
		printf("%d	%f	%f\n", i, x, cx);
	}
}
