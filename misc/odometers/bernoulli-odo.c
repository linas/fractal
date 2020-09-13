/*
 * bernoulli-odo.c
 *
 * Bernoulli odometer
 *
 * Linas Vepstas September 2020
 */

#include <stdbool.h>
#include <math.h>
#include <stdio.h>

/* increment the odometer by one.
 * Unpack the double-precision float into a string of bits,
 * increment, and re-pack into a double.
 */
double increment(double x)
{
	double result = 0.0;
	double bit = 0.5;
	bool done = false;

	for (int i=0; i<56; i++)
	{
		if (done)
		{
			if (0.5 <= x) result += bit;
		}
		else
		{
			if (x < 0.5)
			{
				result += bit;
				done = true;
			}
			else
			{
				/* do nothing here; this rolls over */
			}
		}

		x *= 2.0;
		if (1.0 <= x) x -= 1.0;
		bit *= 0.5;
	}

	return result;
}

int main (int argc, char* argv[])
{
	int nbins=901;

	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);
		double y = increment(x);
		printf("%d	%g	%g\n", i, x, y);
	}
}
