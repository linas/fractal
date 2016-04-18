/*
 * Find zeros of exponential generating function for greatest
 * prime factors.
 *
 * April 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-zerofind.h>
#include <gpf.h>
#include "gpf-gen-bignum.h"

void expo(cpx_t sum, cpx_t z, int nprec)
{
	cpx_gpf_exponential(sum, z, nprec);

	// Divide the returned value by r...
	mpf_t r;
	mpf_init(r);
	cpx_abs(r, z);
	cpx_div_mpf(sum, sum, r);
}

double gpf_bignum_exponential(double r, double theta)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	theta *= 2.0 * M_PI;
	cpx_set_d(z, r*cos(theta), r*sin(theta));

	// cpx_gpf_exponential(sum, z, 20);
	expo(sum, z, 20);

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
	cpx_set_d(e1, 0.25, 0);
	cpx_set_d(e2, 0, 0.25);

	int rc = cpx_find_zero(zero, expo, guess, e1, e2, 20, 50);

	// if rc is not zero, then nothing was found
	if (rc) return;

	double re = cpx_get_re(zero);
	double im = cpx_get_im(zero);

	double r = sqrt(re*re + im*im);
	double t = atan2(im, re) / (2.0*M_PI);

	printf("z = %18.16g * exp(2pi i %18.16g)\n", r, t);
	printf("z = %18.16g + I %18.16g\n", re, im);

	cpx_t check;
	cpx_init(check);
	cpx_gpf_exponential(check, zero, 20);
	double eps_r = cpx_get_re(check);
	double eps_i = cpx_get_im(check);
	double eps = sqrt(eps_r*eps_r + eps_i*eps_i);
	printf("fun at zero = %g\n", eps);
}

void survey(double rmax, double cell_size)
{
	for (double r=cell_size; r<rmax; r += cell_size)
	{
		double prev = 1.0;
		double tprev = 0.0;
		bool reported = false;
		for (double t=0.0; t < 0.5; t += cell_size/r)
		{
			double sample = gpf_bignum_exponential(r, t);

			if (sample < 0.15)
			{
				if (sample < prev)
				{
					prev = sample;
					tprev = t;
					reported = false;
				}
				else if (false == reported)
				{
					printf("---------\n");
					printf("Candidate zero near r=%g t=%g\n", r, tprev);
					fflush(stdout);
					find_zero(r, tprev, cell_size);
					reported = true;
				}
			}
			else
			{
				prev = 1.0;
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
