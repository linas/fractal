
/*
 * Explore sums of binary digits
 *
 * Linas Veptas 2015 August
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sum.h"

int main (int argc, char * argv[])
{
#if 0
	for (double x=0.0; x<1.0; x+= 0.1)
	{
		int xbits[LEN];
		float_to_bitstring(x, xbits);

		for (double y=0.0; y<1.0; y+= 0.08)
		{
			int ybits[LEN];
			float_to_bitstring(y, ybits);

			int bsum[LEN];
			// add_bitstrings(bsum, xbits, ybits);
			reverse_add(bsum, xbits, ybits);
			double sum = bitstring_to_float(bsum);

			double mod = x+y;
			if (1.0 < mod) mod -= 1.0;
			printf("duude %g + %g  = %g \n", x, y, sum);
		}
	}
#endif

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <x>\n", argv[0]);
		exit (-1);
	}
	double x = atof(argv[1]);

	int xbits[LEN];
	float_to_bitstring(x, xbits);

	double ystep = 1.0 / 3197.0;
	for (double y=0.0; y<1.0; y+= ystep)
	{
		int ybits[LEN];
		float_to_bitstring(y, ybits);

		int bsum[LEN];
		// add_bitstrings(bsum, xbits, ybits);
		reverse_add(bsum, xbits, ybits);
		double sum = bitstring_to_float(bsum);

		double mod = x+y;
		if (1.0 < mod) mod -= 1.0;
		printf("%g\t%g\t%g\n", x, y, sum);
	}
}
