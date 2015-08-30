
/*
 * Explore sums of binary digits
 *
 * Linas Veptas 2015 August
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define LEN 60

void float_to_bitstring(double x, int bits[LEN])
{
	x -= floor(x);
	for (int i=0; i<LEN; i++)
	{
		if (0.5 <= x) bits[i] = 1;
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

void add_bitstrings(int sum[LEN], int a[LEN], int b[LEN])
{
	int carry = 0;
	for (int i=LEN-1; i>= 0; i--)
	{
		int s = a[i] + b[i] + carry;
		if (s < 2) { sum[i] = s; carry = 0; }
		else { sum[i] = s-2; carry = 1; }
	}
}

void reverse_add(int sum[LEN], int a[LEN], int b[LEN])
{
	int carry = 0;
	for (int i=0; i<LEN; i++)
	{
		int s = a[i] + b[i] + carry;
		if (s < 2) { sum[i] = s; carry = 0; }
		else { sum[i] = s-2; carry = 1; }
	}
}

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
