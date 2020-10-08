/*
 * extended-debug.C
 *
 * Verify that the extended map is creating valid beta-expansions...
 * Yes, it always was ... but we recorded orbits incorrectly at
 * the branch points. Boo.  But we found and fixed that now.
 * So running this should pass. That is, this is a unit test.
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
		fprintf(stderr, "Usage: %s K maxdepth\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double beta = 2.0*Kay;

	int maxdepth = atof(argv[2]);

	int em = emrun(Kay);
	double a = 0.5;
	double b = a * (1.0 + pow(beta, -em));

#define EPS 1e-16
	int nbits = -log(EPS) / log(beta);

	do_init(nbits);

	mpf_class ebeta = beta;
	double debeta = mpf_get_d(ebeta.get_mpf_t());

	printf("#\n# K=%g em=%d ebeta=%g nbits=%d b=%g\n#\n", Kay, em, debeta, nbits, b);

#define NSAMP 1000
	for (int i=0; i<NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);

		mpf_class ex=x;

		std::vector<std::vector<mpf_class>> orbit_set;
		std::vector<std::vector<bool>> bitset;
		std::vector<std::vector<int>> branch_set;
		std::vector<std::vector<bool>> gamma_set;
		beta_expand(ex, ebeta, em, maxdepth,
		            orbit_set, bitset, branch_set, gamma_set, nbits);

		// printf("initial %g\n", x);

		int ntracks = bitset.size();
		// printf("total of %d tracks\n", ntracks);
		for (int j=0; j<ntracks; j++)
		{
			std::vector<mpf_class> orbit = orbit_set[j];
			std::vector<int> branch_points = branch_set[j];
			std::vector<bool> bits = bitset[j];
			std::vector<bool> gamm = gamma_set[j];

			double y = 0.0;
			double ob = 0.5;
			for (int k=0; k<(int)bits.size(); k++)
			{
				y += bits[k] * ob;
				ob /= beta;
			}
			if (9.0*EPS < fabs(y-x))
			{
				printf("gamma= ");
				print_bits(gamm);
				print_branches(branch_points);
				print_bits(bits);
				printf("j=%d x= %g  expand= %g diff= %g\n", j, x, y, y-x);

				// for (int k=0; k<50; k++)
				//	printf("%d orbit=%g\n", k, mpf_get_d(orbit[k].get_mpf_t()));
				printf("\n");
			}
		}
	}
	printf("Finished %d samples of depth=%d\n", NSAMP, maxdepth);
}

// ================================================================
