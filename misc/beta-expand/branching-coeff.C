/*
 * branching-coeff.C
 * Try to guess what the branching measure is.
 * So far, the attempts to guess all fail, but they come close...
 *
 * Linas Vepstas Oct 2020
 */

#include "sidorov-big.C"

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s y nbits depth\n", argv[0]);
		exit (1);
	}
	double why = atoi(argv[1]);
	int base_nbits = atoi(argv[2]);
	int maxdepth = atoi(argv[3]);
#define MAXDEPTH maxdepth

#define NBINS 503

	printf("#\n# Step function sums at y=%g nbits=%d depth=%d\n#\n",
	      why, base_nbits, maxdepth);

	printf("# column labels\n"
	       "# beta, em, theta(beta/2-y), theta(alpha-y), theta(alpha-beta/2-y),\n"
	       "# integrals for each of the three above.\n#\n");

	for (int i=0; i< NBINS; i++)
	{
		double Kay = (((double) i) + 0.5) / ((double) NBINS);
		Kay = 0.5* (Kay+1.0);
		int em = emrun(Kay);
		double dbeta = 2.0*Kay;
		double dalpha = 0.5 * (1.0 + pow(dbeta, -em));

		int nbits = base_nbits / (4*(2-dbeta)*(dbeta-1));
		do_init(nbits);

		mpf_class beta = dbeta;

		// Orbits of 1.0 aka orbits of "lower bound" beta/2
		// That is, we do the first iteration by hand.
		std::vector<std::vector<mpf_class>> orbit_set;
		std::vector<std::vector<bool>> bitset;
		std::vector<std::vector<int>> branch_set;
		std::vector<std::vector<bool>> gamma_set;
		mpf_class one = dbeta * 0.5;
		beta_expand(one, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, gamma_set, nbits);

		// Orbits of the upper bound, alpha*beta = beta/2 (1+\beta^-m)
		std::vector<std::vector<mpf_class>> erbit_set;
		std::vector<std::vector<bool>> ebitset;
		std::vector<std::vector<int>> ebr_set;
		std::vector<std::vector<bool>> egam_set;

		mpf_class alpha = dbeta * dalpha;
		beta_expand(alpha, beta, em, MAXDEPTH, erbit_set, ebitset, ebr_set, egam_set, nbits);
		int ntracks = bitset.size();
		int etracks = ebitset.size();
		// printf("# ntracks=%d etrack=%d\n", ntracks, etracks);
		if (ntracks < (1<<MAXDEPTH) or etracks < (1<<MAXDEPTH)) continue;

		// printf("# K=%g beta=%g alpha=%g\n", Kay, dbeta, dalpha);
		// printf("# (beta-1)(1-alpha) = %g\n", (dbeta-1.0)*(1.0-dalpha));
		// printf("# 2(beta-1)(1-alpha) = %g\n", 2*(dbeta-1.0)*(1.0-dalpha));

		double y = why;

		double xacc = 0.0;
		double uacc = 0.0;
		double macc = 0.0;

		double Sxacc = 0.0;
		double Suacc = 0.0;
		double Smacc = 0.0;

		double bpn = 1.0;

		// Get a reasonable amount of precision for small beta.
		int kmax=-log(1e-6)/log(dbeta);
		if ((int) bitset[0].size() < kmax) kmax=bitset[0].size() -1;

		// Sum over powers 1/beta^k
		for (int k=0; k<kmax; k++)
		{
			// Average over paths gamma
			for (int j=0; j<ntracks; j++)
			{
				// std::vector<bool>& bits = bitset[j];
				std::vector<mpf_class>& orbit = orbit_set[j];
				std::vector<int>& branch_points = branch_set[j];

				size_t nb = branch_points.size();
				size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
				if (orbit.size() <= norb) norb = orbit.size() -1;
				if ((int) norb < k) continue;

				// std::vector<bool>& ebits = ebitset[j];
				std::vector<mpf_class>& erbit = erbit_set[j];
				// std::vector<int>& ebranches = branch_set[j];

				double x = mpf_get_d(orbit[k].get_mpf_t());
				double xu = mpf_get_d(erbit[k].get_mpf_t());

				if (y < x) xacc += bpn;
				if (y < xu) uacc += bpn;
				if (y < xu - x) macc += bpn;

				// These are the integrals - the normalization.
				Sxacc += x*bpn;
				Suacc += xu*bpn;
				if (0 < xu-x) Smacc += (xu-x)*bpn;
			}
			bpn /= dbeta;
		}
		xacc /= ntracks;
		uacc /= ntracks;
		macc /= ntracks;
		Sxacc /= ntracks;
		Suacc /= ntracks;
		Smacc /= ntracks;
		printf("%g	%d	%g	%g	%g	%g	%g	%g\n",
		       dbeta, em, xacc, uacc, macc, Sxacc, Suacc, Smacc);
	}
}

// ================================================================
