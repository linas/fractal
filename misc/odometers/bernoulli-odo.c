/*
 * bernoulli-odo.c
 *
 * Bernoulli odometer
 *
 * Linas Vepstas September 2020
 */

#include <complex.h>
#include <math.h>

/* increment the odometer by one.
 * Unpack the double-precision float into a string of bits,
 * increment, and re-pack into a double.
 */
double increment(double x)
{
	double result;

	for (int i=0; i<56; i++)
	{
		bool bit = (x

		x *= 2.0
		if (1.0 <= x) x -= 1.0;
	}

	return result;
}
