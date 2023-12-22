
/*
 * Simple graphs of the bit-wise operations.
 *
 * Jan 2018 Linas Vepstas
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bitops.h"

int clamp(int n)
{
	return 2*(n/2);

	return n/2;

	return n%2;

	int max = 2;
	if (max < n) return max;
	return n;
}

int main(int argc, char* argv[])
{
	double K1 = atof(argv[1]);
	double K2 = atof(argv[2]);
	double K3 = atof(argv[3]);
	double K4 = atof(argv[4]);

#define NBINS 979
	for (int i=0; i< NBINS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NBINS);
/*
		double k1x = add_xor(K1, x);
		double k2x = add_xor(K1, k1x);
		double k3x = add_xor(K3, x);
		double k4x = add_xor(K4, x);
*/
/*
		K3 = add_xor(K1,K2);
		double k1x = mult_xor(K1, x);
		double k2x = mult_xor(K2, x);

		double k3x = mult_xor(K3, x);
		double k4x = mult_xor(K4, x);
		k4x = add_xor(k1x, k2x);
*/

/*
		double k1x = mult_xor(K1, x);
		double k2x = mult_xor(K1, k1x);
		double k3x = mult_xor(K1, k2x);
		double k4x = mult_xor(K1, k3x);
*/
		double k1x = mangle_multiply(K1, x, clamp);
		double k2x = mangle_multiply(K2, x, clamp);
		double k3x = mangle_multiply(K3, x, clamp);
		double k4x = mangle_multiply(K4, x, clamp);
k3x = mult_xor(K1, x);
k4x = mult_xor(K2, x);

		printf ("%d	%g	%g	%g	%g	%g\n", i, x, k1x, k2x, k3x, k4x);
	}
}
