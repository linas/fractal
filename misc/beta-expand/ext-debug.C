/*
 * ext-debug.C
 * Verify that the extended map is creating valid beta-expansions...
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "emrun.C"

std::vector<bool> bits;

void prt_bits()
{
	printf("len=%lu ", bits.size());
	for (size_t j=0; j<bits.size(); j++)
		printf("%d", (int) bits[j]);
	printf("\n");
}

double tee(double x, double beta)
{
	if (x<0.5) return beta*x;
	if (x<1.0) return beta*(x-0.5);
	return beta*(x-1.0);
}

double maybe(double x, double beta)
{
	long int r = random();
	if (r < RAND_MAX/2)
	{
		bits.push_back(0);
		return beta*x;
	}
	bits.push_back(1);
	return tee(x, beta);
}

double tau(double x, double beta, double a, double b)
{
	if (a < x and x < b) return maybe(x, beta);
	if (a+0.5 < x and x < b+0.5) return maybe(x, beta);
	if (a+1.0 < x and x < b+1.0) return maybe(x, beta);

	if (x<0.5) { bits.push_back(0); return beta*x; }
	if (x<1.0) { bits.push_back(1); return beta*(x-0.5); }

	// fprintf(stderr, "fail x=%g\n", x);
	// exit(1);
	return beta*(x-1.0);
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double beta = 2.0*Kay;

	int em = emrun(Kay);
	double a = 0.5/beta;
	double b = a * (1.0 + pow(beta, -em));
	printf("#\n# K=%g m=%d\n#\n", Kay, em);


#define NSAMP 8
	for (int i=0; i<NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);

		bits.clear();
		for (int j=0; j<50; j++)
		{
			x = tau(x, beta, a, b);
		}

		prt_bits();

		double y = 0.0;
		double ob = 1.0 / beta;
		for (int j=0; j<50; j++)
		{
			y += bits[j] * ob;
			ob /= beta;
		}
		printf("%g	%g\n", x, y);
	}
}

// ================================================================
