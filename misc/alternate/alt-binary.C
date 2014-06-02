
/**
 * Alternating sums inspired by differential oparators
 *
 * Linas Vepstas
 * June 2014
 */

#include <math.h>

#define BITLEN 56
typedef char[BITLEN] bitstring;

bitstring long to_bits(double x)
{
	bitstring bits;
	x -= floor(x);

	for (int i=0; i<BITLEN; i++)
	{
		if (0.5 < x) bits[i] = 1;
		else bits[i] = 0;
	}

	return bits;
}
