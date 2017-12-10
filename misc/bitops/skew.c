
/*
 * skew.c
 *
 * Skew takago map
 * Dec 2017
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double skew_triangle(double x, double K)
{
	x -= floor(x);
	if (x < 0.5*K) return 2.0 * x / K;
	return (1.0 - x) / (1.0 - 0.5*K);
}

double skew_takagi(double x, double K, double w, int cnt)
{
	if (cnt < 0)
		return skew_triangle(x,K);

	if (x < 0.5*K) 
	{
		return 2.0 * x / K + w*skew_takagi(x/(0.5*K), K, w, cnt-1);
	}
	else
	{
		return (1.0 - x) / (1.0 - 0.5*K) + w*skew_takagi(x/(1-0.5*K), K, w, cnt-1);
	}
}

int main (int argc, char* argv[])
{
	double K = atof(argv[1]);
	double w = atof(argv[2]);

#define NPTS 803
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double y = skew_takagi(x, K, w, 2);
		printf("%d	%g %g\n", i, x, y);
	}
}
