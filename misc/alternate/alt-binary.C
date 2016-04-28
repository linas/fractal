
/**
 * Alternating sums inspired by differential operators.
 * Viz. just like in homology: the boundary operator drops the n'th
 * index, and sums, alternating the sign as (-1)^n.  We do the same
 * here, where the index is the bit position in a binary string.
 * The result is a fractal .. which I've seen before, but where? why?
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

	for (double x=0.0; x<1.0; x+=0.00212345678)
	{
		double y = alternate(x);

		// z equals y exactly ...
		// double z = 4.0*alternate(0.25*x);

		// W equals y exactly.
		// double w = (x>0.5)? -0.5 : 0.5;
		// w += alternate(x+0.5);

		double u = alternate(y);

		printf("%16.14g	%16.14g	%16.14g\n", x, y, u);
	}
}
