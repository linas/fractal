/*
 * sidorov-measure.C
 * Take a stab at the a-priori measure.
 *
 * Linas Vepstas Sept 2020
 */

#include "sidorov-big.C"

int main (int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, "Usage: %s K nbits depth\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	int nbits = atoi(argv[2]);
	int maxdepth = atoi(argv[3]);

	do_init(nbits);

	mpf_class beta;
	// make_random_bitsequence(beta, 2.0*Kay, nbits, NBINS);
	beta = 2.0*Kay;

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d nbits=%d\n#\n", Kay, em, nbits);

#define MAXDEPTH maxdepth
	mpf_class ex = 1;

	std::vector<std::vector<mpf_class>> orbit_set;
	std::vector<std::vector<bool>> bitset;
	std::vector<std::vector<int>> branch_set;
	beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

	int ntracks = bitset.size();
	for (int j=0; j<ntracks; j++)
	{
		std::vector<mpf_class> orbit = orbit_set[j];
		std::vector<int> branch_points = branch_set[j];

		size_t nb = branch_points.size();
		size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
		if (orbit.size() <= norb) norb = orbit.size() -1;
		for (size_t k=1; k<=norb; k++)
		{
			double x = mpf_get_d(orbit[k].get_mpf_t());
		}
	}
}

// ================================================================
