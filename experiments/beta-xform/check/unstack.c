/*
 * unstack.c
 *
 * Verify to unwrapped recursion relations for generalized
 * stretch-cut-stack map. The are the ones after noticing that
 * the unwrap simplifies even more.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ==============================================================

// Arbitrary function
double nu(double x)
{
	return 1.0;
	// return x-0.5;

	// Bernoulli poly B_2
	// return x*x - x  + 1.0 / 6.0;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// ==============================================================

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
int b_k(double beta, int k)
{
	double tk = t_k(beta, k);
	if (beta*tk >= 1.0) return 1;
	return 0;
}

// Return the q_k constant
double q_k(double beta, int k)
{
	double sum = 0.0;
	double bej = 1.0;
	for (int j=0; j<k; j++)
	{
		sum += b_k(beta, j) / bej;
		bej *= beta;
	}
	return sum * bej;
}

// ==============================================================

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);

	for (int k=1; k<= n; k++)
	{
		double qk = q_k(beta, k);
		double bkp1 = pow(beta, k+1);
		double tk = t_k(beta, k+1);
		printf("%d	%g	%g	%g\n", k, qk, bkp1-qk, tk);
	}
}
