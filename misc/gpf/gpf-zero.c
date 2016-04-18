/*
 * Find zeros of exponential generating function for greatest
 * prime factors.
 *
 * April 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gpf.h>
#include "gpf-gen-bignum.h"

double gpf_bignum_exponential(double r, double theta)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	theta *= 2.0 * M_PI;
	cpx_set_d(z, r*cos(theta), r*sin(theta));

	cpx_gpf_exponential(sum, z, 20);

	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);
	return rv;
}

void survey(double rmax)
{
	for (double r=1.0; r<rmax; r += 1.0)
	{
		for (double t=0.0; t < 1.0; t += 1.0/r)
		{
			double sample = gpf_bignum_exponential(r, t);
			sample /= r;

			if (sample < 1.0)
			{
				printf("duude got one here: %g %g\n", r, t);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <rmax>\n", argv[0]);
		exit(1);
	}

	double rmax = atof(argv[1]);
	survey(rmax);
	fflush(stdout);
}
