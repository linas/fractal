/*
 * midmap.c
 *
 * Plot the left-right moves of the midpoint.
 * As the midpoint is iterated, it makes a sequence of left-right
 * moves. Convert these into a binary number, and graph it.
 * It looks like some deformed compressor function, not particularly
 * enlightening.
 *
 * February 2018
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double midpoint_orbit(double K)
{
	double map = 0.0;
	double hn = 0.5;

	double m = K;
	for (int i=0; i< 50; i++)
	{
		if (m <= 0.5) m = 2.0*K*m;
		else
		{
			m = 2.0*K*m - K;
			map += hn;
		}
		hn *= 0.5;
	}
	return map;
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s npts\n", argv[0]);
		exit(1);
	}
	int npts = atoi(argv[1]);

	for (int i=0; i< npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		double K = 0.5 + 0.5*x;
		double y = midpoint_orbit(K);
		y = 2.0*(y-0.5);
		printf("%d	%g	%g\n", i, 2.0*K, y);
	}
}
