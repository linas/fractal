
/*
 * ftag.c
 *
 * Rebuild the Farey dist as a Cantor set
 * This time, as a log of Takagi thing.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* A stupa is a step pyramid */
double stupa (double x)
{
	int n;

	x -= floor(x);
	if (0.5 < x) x = 1.0-x;
	if (0.0 == x) return 0.0;
	
	n = 0;
	while(1)
	{
		if (0.25 < x) break;
		n++;
		x *= 2.0;
	}
	return n;
}

double taga(double x, double weight)
{
	int k;

	double acc = 0.0;
	double tk = 1.0;
	double tw = 1.0;
	for (k=0; k<20; k++)
	{
		tak
	}

}
