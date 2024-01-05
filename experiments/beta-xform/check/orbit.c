/*
 * orbit.c
 *
 * Manual exploration of individual orbits.
 * Jan 2024
 */

#include <stdio.h>
#include <stdlib.h>

void orbit(double beta)
{
	double mid = 0.5* beta;
	for (int i=0; i<6; i++)
	{
		printf("%d	%g\n", i, mid);
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);

	orbit(beta);
}
