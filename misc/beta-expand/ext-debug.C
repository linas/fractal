/*
 * ext-debug.C
 *
 * Verify that the extended map is creating valid beta-expansions...
 * and now that we bug-fixed everything ... it is! Hooray!
 * This is a unit test.  It should pass.
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "emrun.C"

std::vector<bool> bits;
std::vector<bool> gamm;
int lvl = 0;

void prt_bits(std::vector<bool> stuff)
{
	printf("len=%lu ", stuff.size());
	for (size_t j=0; j<stuff.size(); j++)
		printf("%d", (int) stuff[j]);
	printf("\n");
}

double tee(double x, double beta)
{
	if (x<0.5) return beta*x;
	return beta*(x-0.5);
	// if (x<1.0) return beta*(x-0.5);
	// return beta*(x-1.0);
}

double maybe(double x, double beta)
{
	long int r = random();
	if (r < RAND_MAX/2)
	{
		// printf(" went low\n");
		gamm.push_back(0);
		bits.push_back(0);
		return beta*x;
	}

	gamm.push_back(1);
	// printf(" went high\n");
	if (x<0.5) { bits.push_back(0); return beta*x; }
	if (x<1.0) { bits.push_back(1); return beta*(x-0.5); }

	fprintf(stderr, "------------fool x=%g\n", x);
	exit(1);
	bits.push_back(1);
	return beta*(x-1.0);
}

double tau(double x, double beta, double a, double b)
{
	// printf("x= %g cyl=%d lvl=%d\n", x, a < x and x < b, lvl);
	lvl++;
	if (a < x and x < b)
	{
		return maybe(x, beta);
	}

	if (x<0.5) { bits.push_back(0); return beta*x; }
	bits.push_back(1);
	return beta*(x-0.5);
}

// Given branch points gamma, recreate tau.
// This is the one-liner formula in the paper.
// It works.
double tau_binary(double x, double beta, double a, double b,
                  std::vector<bool>& gammy)
{
	double y = tee(x, beta);
	if (a < x and x < b)
	{
		bool g = gammy.front();
		gammy.erase(gammy.begin());
		if (not g) y += 0.5*beta;
		// printf("bing %d %lu\n", g, gammy.size());
	}
	return y;
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
	double a = 0.5;
	double b = a * (1.0 + pow(beta, -em));
	printf("#\n# K=%g m=%d\n#\n", Kay, em);

#define EPS 1e-14
	int nbits = -log(EPS) / log(beta);

#define NSAMP 1000000
	for (int i=0; i<NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);

		// Record an orbit
		bits.clear();
		gamm.clear();
		double z = x;
		std::vector<double> orbit;
		for (int j=0; j<nbits; j++)
		{
			orbit.push_back(z);
			z = tau(z, beta, a, b);
		}

		// Verify alternate orbit ... the single-line reconstructed orbit.
		z = x;
		std::vector<bool> gammy = gamm;
		for (int j=0; j<nbits; j++)
		{
			if (EPS < fabs(z-orbit[j]))
			{
				// Rounding errors will accumulate ... ignore those.
				if (10*EPS < fabs(z-orbit[j]))
					printf("failed reconstruction at %d %g %g %g\n",
					        j, z, orbit[j], z-orbit[j]);
				break;
			}
			z = tau_binary(z, beta, a, b, gammy);
		}
		// printf("\n");

		double y = 0.0;
		double ob = 0.5;
		for (int j=0; j<nbits; j++)
		{
			y += bits[j] * ob;
			ob /= beta;
		}
		if (9.0*EPS < fabs(y-x))
		{
			printf("gamma= ");
			prt_bits(gamm);
			prt_bits(bits);
			printf("i=%d x= %g  expand= %g diff= %g\n", i, x, y, y-x);
			printf("\n");
		}
	}

	printf("Complete running %d samples\n", NSAMP);
}

// ================================================================
