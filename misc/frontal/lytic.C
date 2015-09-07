/*
 *
 * Numerical exploration of the analytic continuation
 * given in frontal.lyx
 *
 * September 2015
 */

#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
// #include <complex>

// using std::complex;

#include "lytic.h"

/*
 * Return the count c_k of the k'th consequtive sequence of bits.
 * Sometimes called the "continued fraction expansion" or "Denjoy sum".
 *
 * Note the normalization here: c_1 is the number of leding zero's in
 * the expansion. It can be zero.
 */
int bitcount(int k, int bits[LEN])
{
	int bp = 0;
	int parity = 0;
	while (k != 1)
	{
		while (parity == bits[bp] and bp < LEN)
		{
			bp++;
		} 
		k -= 1;
		parity = 1 - parity;
	}
	
	int count = 0;
	while (parity == bits[bp] and bp < LEN)
	{
		count++;
		bp++;
	} 

	return count;
}

int xbitcount(int k, double x)
{
	int bits[LEN];
	float_to_bitstring(x, bits);
	return bitcount(k, bits);
}

double complex count_extend(int k, int c_k, int c_1, double complex u)
{
	if (2 < k and k%2 == 0)
		return c_k * u / (1.0 - u);
	if (2 < k and k%2 == 1)
		return c_k * (1.0 - 2.0 * u) / (1.0 - u);

	if (2 == k and 0 < c_1)
		return c_k * u / (1.0 - u);
	if (2 == k and 0 == c_1)
		return (u * (1.0 - 2.0 * u) / (1.0 - u)) + (c_k - 1) * u / (1.0 - u);

	if (1 == k and 0 < c_1)
		return (u * (1.0 - 2.0 * u + 2.0 *u *u) / (1.0 - u))
			+ (c_k - 1) * (1.0 - 2.0 * u) / (1.0 - u);
	return 0.0;
}

double complex count_extend(int k, int bits[LEN], double complex u)
{
	int c_1 = bitcount(1, bits);
	int c_k = bitcount(k, bits);
	return count_extend(k, c_k, c_1, u);
}

double complex count_extend(int k, double x, double complex u)
{
	int bits[LEN];
	float_to_bitstring(x, bits);
	return count_extend(k, bits, u);
}

