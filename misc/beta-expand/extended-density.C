/*
 * extended-density.C
 * Compute the extened density measure
 *
 * Linas Vepstas Sept 2020
 */

#include "sidorov-big.C"

void extended_measure(double dbeta,
                      double* array,  // where to store the density
                      int NBINS,      // length of above
                      int nbits)      // precision
{
	// Adjust precision for beta; this is the U-shaped curve.
	nbits *= 1.0/(4.0*(dbeta-1.0)*(2.0-dbeta));
	if (1.6 < dbeta) nbits *= 0.4/(2.0-dbeta);  // approx fit to high range
#define MAXBITS 1200
	if (nbits < 0) nbits=MAXBITS;
	if (MAXBITS < nbits) nbits=MAXBITS;
	do_init(nbits);

	mpf_class beta;
	make_random_bitsequence(beta, dbeta, nbits, NBINS);
	int em = emrun(0.5*dbeta);

	/* clear out the row */
	for (int j=0; j<NBINS; j++) array[j] = 0.0;

	fprintf(stderr, "Working K=%g beta=%g nbits=%d\n", dbeta/2.0, dbeta, nbits);

	for (int ibin=0; ibin<NBINS; ibin++)
	{
		// if (ibin%100 ==0) fprintf(stderr, "# orbits done %d of %d\n", ibin, NBINS);
		// fprintf(stderr, "# orbits done %d of %d\n", ibin, NBINS);
		double x = ((double) ibin)/ ((double) NBINS);

#define NSAMP 32
		for (int nsamp=0; nsamp<NSAMP; nsamp++)
		{
			mpf_class ex;
			make_random_bitsequence(ex, x, nbits, NBINS);

#define MAXDEPTH 6
			std::vector<std::vector<mpf_class>> orbit_set;
			std::vector<std::vector<bool>> bitset;
			std::vector<std::vector<int>> branch_set;
			beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

			// Compute a histogram of the orbits. But do it
			// by summing only up to the last branch-point.
			// (else the greedy expansion will dominate).
			int ntracks = bitset.size();
			for (int j=0; j<ntracks; j++)
			{
				std::vector<mpf_class> orbit = orbit_set[j];
				std::vector<int> branch_points = branch_set[j];
#define PHI (0.5 * (sqrt(5.0) + 1.0))
#define SCALE (0.5 * (PHI + 1.0))
				size_t nb = branch_points.size();
				size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
				if (orbit.size() <= norb) norb = orbit.size() -1;
				for (size_t k=1; k<=norb; k++)
				{
					double x = mpf_get_d(orbit[k].get_mpf_t());
					if (10.0 < x) continue;  // wtf!???
					int bin = x * NBINS / SCALE;
					if (NBINS <= bin) bin=NBINS-1;
					array[bin] += 1.0;
				}
			}
		}
	}

	// Normalize
	double nobs = 0;
	for (int j=0; j<NBINS; j++)
		nobs += array[j];
	for (int j=0; j<NBINS; j++)
		array[j] *= NBINS/nobs;
}

// ================================================================
