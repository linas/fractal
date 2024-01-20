/*
 * unrolling.c
 * Verify to unrolled recursion relations for generalized
 * stretch-cut-stack map.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Return endpoint iterate.
double t_k(double beta, int k)
{
	double tk = 1.0;
	for (int i=0; i<k; i++)
	{
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
	}
	return tk;
}

// Return midpoint iterate digit b_k = d_k(1/2)=theta(beta t_k-1)
double b_k(double beta, int k)
{
	double tk = t_k(beta, k);
	if (beta*tk >= 1.0) return 1.0;
	return 0.0;
}

// Arbitrary function
double nu(double x)
{
	return 1.0;
}

// Return the g_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double g_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (0 == k) return nu(x);

	if (1 == k)
	{
		printf("not done yet\n");
		return 1.0;
	}
	double bkm1 = b_k(beta, k-1);
	double arg = (x + bkm1)/ beta;
	double guh = g_n_k(beta, n-1, k-1, arg);
	return guh / beta;
}

// Return the g_n_1 constant from the "generalized stretch-cut-stack"
// section of paper.
double g_n_1(double beta, int n)
{
	return 0.0;
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta \n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);

	g_n_k(beta, 3, 3, 0.3);
}
