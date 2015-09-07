/*
 *
 * Numerical exploration of the analytic continuation
 * given in frontal.lyx
 *
 * September 2015
 */

#include <stdio.h>
#include <stdlib.h>

#include "../sum/sum.h"

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

int main(int argc, char *argv[])
{
	double x = 0.49999999999;

	x = atof(argv[1]);

	for (int k=1; k<10; k++)
	{
		int c_k = xbitcount(k, x);
		printf("duuude its %d %d\n", k, c_k);
	}
}
