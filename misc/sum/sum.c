
/*
 * Explore sums of binary digits
 *
 * Linas Veptas 2015 August
 */

#include <math.h>
#include <stdio.h>

#define LEN 60

void float_to_bitstring(double x, int bits[LEN])
{
	x -= floor(x);
	for (int i=0; i<LEN; i++)
	{
		if (0.5 < x) bits[i] = 1;
		else bits[i] = 0;
		x *= 2.0;
		x -= floor(x);
	}
}

double bitstring_to_float(int bits[LEN])
{
	double x = 0.0;
	double tn = 0.5;
	for (int i=0; i<LEN; i++)
	{
		if (1 == bits[i]) x += tn;
		tn *= 0.5;
	}
	return x;
}

int main (int argc, char * argv[])
{
	for (double x=0.0; x<1.0; x+= 0.0123)
	{
		int bits[LEN];
		float_to_bitstring(x, bits);
		double y = bitstring_to_float(bits);

		printf("duude %g %g \n", x, x-y);
	}

}
