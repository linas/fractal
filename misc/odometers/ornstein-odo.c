/*
 * ornstein-odo.c
 *
 * Ornstein odometer
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
double orn_increment(double x)
{
	double result = 0.0;
	double fact = 0.5;
	bool done = false;

	// factorial(19) is just a little bigger than 2^56
	for (int i=0; i<19; i++)
	{
		double bit = floor (x * (i+2));
		if (done)
		{
			result += bit * fact;
		}
		else
		{
			if (bit < i+1)
			{
				result += (bit + 1.0) * fact;
				done = true;
			}
			else
			{
				/* do nothing here; this rolls over */
			}
		}

		x *= (double) i+2;
		x = floor(x);
		fact /= i+3;
	}

	return result;
}

int main (int argc, char* argv[])
{
	int nbins=901;

	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);
		double y = orn_increment(x);
		double y2 = orn_increment(y);
		double y3 = orn_increment(y2);
		double y4 = orn_increment(y3);
		printf("%d	%g	%g	%g	%g	%g\n", i, x, y, y2, y3, y4);
	}
}
