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
	printf("#\n# K=%g m=%d b=%g\n#\n", Kay, em, b);

#define EPS 1e-20
	int nbits = -log(EPS) / log(beta);

	do_init(nbits);

	mpf_class ebeta;
	make_random_bitsequence(ebeta, beta, nbits, 1);


#define NSAMP 1
	for (int i=0; i<NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);

		mpf_class ex;
		make_random_bitsequence(ex, x, nbits, 1);

		std::vector<std::vector<mpf_class>> orbit_set;
		std::vector<std::vector<bool>> bitset;
		std::vector<std::vector<int>> branch_set;
		beta_expand(ex, ebeta, em, maxdepth, orbit_set, bitset, branch_set, nbits);


		int ntracks = bitset.size();
		printf("totoal of %d tracks\n", ntracks);
		for (int j=0; j<ntracks; j++)
		{
			std::vector<mpf_class> orbit = orbit_set[j];
			std::vector<int> branch_points = branch_set[j];
			std::vector<bool> bits = bitset[j];

			size_t nb = branch_points.size();
			size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
			if (orbit.size() <= norb) norb = orbit.size() -1;
			printf("orbit=%lu\n", norb);

			double y = 0.0;
			double ob = 0.5;
			for (int k=0; k<nbits; j++)
			{
				y += bits[k] * ob;
				ob /= beta;
			}
			// if (9.0*EPS < fabs(y-x))
			{
				// prt_bits();
				printf("i=%d x= %g  expand= %g diff= %g\n", i, x, y, y-x);
				printf("\n");
			}
		}
	}
}

// ================================================================
