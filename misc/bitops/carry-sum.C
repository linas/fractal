/*
 * carry-sum.C
 *
 * Carry-bit dynamics, for sums.
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

void to_bits(char* bits, double x, int len)
{
	// convert x into bit-string.
	for (int i=0; i<len; i++)
	{
		if (0.5 <= x)
		{
			x -= 0.5;
			bits[i] = 1;
		}
		else bits[i] = 0;
		x *= 2;
	}
}

double to_float(char* bits, int len)
{
	// convert bit-string to float;
	double x = 0.0;
	double tn = 0.5;
	for (int i=0; i<len; i++)
	{
		if (bits[i]) x += tn;
		tn /= 2;
	}
	return x;
}

double add_carry(double x, double y)
{
	// convert x and y into bit-strings.
	char xbits[50];
	char ybits[50];
	to_bits(xbits, x, 50);
	to_bits(ybits, y, 50);

	char cbits[50];
	for (int i=50; 0<i; i--)
	{
		int sum = xbits[i]+ybits[i]+cbits[i];
		if (1 < sum) cbits[i-1] = 1;
		else cbits [i-1] = 0;
	}

	// Reconstruct x in a mashed bernoulli sequence.
	return to_float(cbits, 50);
}

static void mapping_diagram (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double y,
                             int itermax,
                             double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double deltax = x_width / array_size;
	double startx = x_center - 0.5 * x_width;
	double x = startx;
	for (int j=0; j<array_size; j++)
	{
		array[j] = add_carry(x,y);
		x += deltax;
	}
}

DECL_MAKE_BIFUR(mapping_diagram)
