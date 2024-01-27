/*
 * unutil.c
 *
 * Utilities for the un-series.
 *
 * January 2024
 */

#ifndef __UNUTIL_C__
#define __UNUTIL_C__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

// Return the Golfond-Parry normalization constant capital F.
double fnorm(double beta)
{
	double sum = 0.0;
	double bei = 1.0;
	for (int i=0; i<2000; i++)
	{
		sum += t_k(beta, i) / bei;
		bei *= beta;
		if (1.0 < bei * 1e-17) break;
	}
	return sum;
}

#endif // __UNUTIL_C__
