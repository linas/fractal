/*
 * Explore sums of binary digits
 *
 * Linas Veptas 2015 August
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sum.h"

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

double rev_add(double x, double y)
{
	int xbits[LEN];
	float_to_bitstring(x, xbits);
	int ybits[LEN];
	float_to_bitstring(y, ybits);

	int bsum[LEN];
	reverse_add(bsum, xbits, ybits);

	double sum = bitstring_to_float(bsum);

	return sum;
}
