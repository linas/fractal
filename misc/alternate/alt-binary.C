
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

void drop_bit(const bitstring &from, bitstring &to, int place)
{
	for (int i=0, j=0; i<BITLEN; i++, j++)
	{
		if (i == place) i++;
		to[j] = from[i];
	}
}

double alternate(double x)
{
	bitstring bs;
	to_bits(bs, x);
	double acc = 0.0;
	double sgn = 1.0;

	for (int i=0; i<BITLEN; i++)
	{
		bitstring drop;
		drop_bit(bs, drop, i);
		double y = from_bits(drop);
		acc += sgn * y;
		sgn = - sgn;
	}

	return acc;
}

main(int argc, char* argv[])
{

	for (double x=0.0; x<1.0; x+=0.00312345678)
	{
		double y = alternate(x);
		printf("%16.14g	%16.14g\n", x, y);
	}
}


