/*
 * Find zeros of exponential generating function for greatest
 * prime factors.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-zerofind.h>
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

void find_zero(double rguess, double tguess, double cell_size)
{
	cpx_t zero;
	cpx_init(zero);

	cpx_t guess;
	cpx_init(guess);
	tguess *= 2.0 * M_PI;
	cpx_set_d(guess, rguess*cos(tguess), rguess*sin(tguess));

	cpx_t e1, e2;
	cpx_init(e1);
	cpx_init(e2);
	cpx_set_d(e1, cell_size, 0);
	cpx_set_d(e2, 0, cell_size);

	cpx_find_zero(zero, cpx_gpf_exponential, guess, e1, e2, 20, 50);

	double re = cpx_get_re(zero);
	double im = cpx_get_im(zero);

	printf("found one at %18.16g %18.16g\n", re, im);
}

void survey(double rmax, double cell_size)
{
	for (double r=cell_size; r<rmax; r += cell_size)
	{
		for (double t=0.0; t < 0.5; t += cell_size/r)
		{
			double sample = gpf_bignum_exponential(r, t);
			sample /= r;

			if (sample < 0.5)
			{
				printf("Candidate zero near r=%g t=%g\n", r, t);
				fflush(stdout);
				find_zero(r, t, cell_size);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <rmax> <cell-size>\n", argv[0]);
		exit(1);
	}

	double rmax = atof(argv[1]);
	double cell_size = atof(argv[2]);
	survey(rmax, cell_size);
}
