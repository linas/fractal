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

	std::vector<std::vector<mpf_class>> orbit_set;
	std::vector<std::vector<bool>> bitset;
	std::vector<std::vector<int>> branch_set;
	mpf_class one = 1;
	beta_expand(one, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

	int ntracks = bitset.size();

#define NBINS 402
	double meas[NBINS];

	for (int i=0; i< NBINS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NBINS);
		double acc = 0.0;

	for (int k=0; k<20; k++)
	{
		for (int j=0; j<ntracks; j++)
		{
			std::vector<mpf_class> orbit = orbit_set[j];
			std::vector<int> branch_points = branch_set[j];

			size_t nb = branch_points.size();
			size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
			if (orbit.size() <= norb) norb = orbit.size() -1;
			if ((int) norb < k) continue;

			double x = mpf_get_d(orbit[k].get_mpf_t());
		}
	}
}

// ================================================================
