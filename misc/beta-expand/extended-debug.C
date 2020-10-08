/*
 * ext-debug.C
 *
 * Verify that the extended map is creating valid beta-expansions...
 * and now that we bug-fixed everything ... it is! Hooray!
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "sidorov-big.C"


int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double beta = 2.0*Kay;

	int maxdepth=5;

	int em = emrun(Kay);
	double a = 0.5;
	double b = a * (1.0 + pow(beta, -em));

#define EPS 1e-20
	int nbits = -log(EPS) / log(beta);

	do_init(nbits);

	mpf_class ebeta = beta;
	double debeta = mpf_get_d(ebeta.get_mpf_t());

	printf("#\n# K=%g em=%d ebeta=%g nbits=%d b=%g\n#\n", Kay, em, debeta, nbits, b);

#define NSAMP 1
	for (int i=0; i<NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);

		mpf_class ex=x;

		std::vector<std::vector<mpf_class>> orbit_set;
		std::vector<std::vector<bool>> bitset;
		std::vector<std::vector<int>> branch_set;
		beta_expand(ex, ebeta, em, maxdepth, orbit_set, bitset, branch_set, nbits);

		printf("initial %g ex=%g\n", x, mpf_get_d(ex.get_mpf_t()));

		int ntracks = bitset.size();
		printf("totoal of %d tracks\n", ntracks);
		for (int j=0; j<ntracks; j++)
		{
			std::vector<mpf_class> orbit = orbit_set[j];
			std::vector<int> branch_points = branch_set[j];
			std::vector<bool> bits = bitset[j];

			double y = 0.0;
			double ob = 0.5;
			for (int k=0; k<(int)bits.size(); k++)
			{
				y += bits[k] * ob;
				ob /= beta;
			}
			// if (9.0*EPS < fabs(y-x))
			{
				// prt_bits();
				printf("j=%d x= %g  expand= %g diff= %g\n", j, x, y, y-x);
				printf("\n");
			}
		}
	}
}

// ================================================================
