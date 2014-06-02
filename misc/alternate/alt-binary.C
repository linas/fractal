
/**
 * Alternating sums inspired by differential oparators
 *
 * Linas Vepstas
 * June 2014
 */

#include <math.h>
#include <stdio.h>

#define BITLEN 56
typedef char bitstring[BITLEN];

void to_bits(bitstring bits, double x)
{
	x -= floor(x);

	for (int i=0; i<BITLEN; i++)
	{
		if (0.5 < x) { bits[i] = 1; x -= 0.5; }
		else bits[i] = 0;
		x *= 2.0;
	}
}

double from_bits(bitstring bits)
{
	double val = 0.0;
	double bv = 0.5;
	for (int i=0; i<BITLEN; i++)
	{
		if (bits[i]) val += bv;
		bv *= 0.5;
	}
	return val;
}

main(int argc, char* argv[])
{

	double x = 0.1233456789;
	bitstring bs;
	to_bits(bs, x);
	double y = from_bits(bs);

	printf("duue its %g\n", y);
}


