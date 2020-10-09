/*
 * sidorov-measure.C
 * Try to guess what the extended measure is.
 * So far, the attempts to guess all fail, but they come close...
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

	double dbeta = 2.0*Kay;

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d nbits=%d\n#\n", Kay, em, nbits);

#define MAXDEPTH maxdepth

	// Orbits of 1.0 aka orbits of beta/2
	// That is, we do the first iteration by hand.
	std::vector<std::vector<mpf_class>> orbit_set;
	std::vector<std::vector<bool>> bitset;
	std::vector<std::vector<int>> branch_set;
	std::vector<std::vector<bool>> gamma_set;
	mpf_class one = Kay;
	mpf_class beta = dbeta;
	beta_expand(one, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, gamma_set, nbits);

	// Orbits of the upper bound, 0.5 (1+\beta^-m)
	// Do the first iteration by hand.
	std::vector<std::vector<mpf_class>> erbit_set;
	std::vector<std::vector<bool>> ebitset;
	std::vector<std::vector<int>> ebr_set;
	std::vector<std::vector<bool>> egam_set;

	// XXX FIXME I'm confused why is MAXDEPTH off by one?
	mpf_class edge = Kay * pow(dbeta, -em);
	beta_expand(edge, beta, em, MAXDEPTH+1, erbit_set, ebitset, ebr_set, egam_set, nbits);

	int ntracks = bitset.size();
	int etracks = ebitset.size();
	printf("# ntracks=%d etrack=%d\n", ntracks, etracks);

	std::vector<mpf_class> parry_orbit;
	std::vector<bool> parry_bits;
	beta_sequence(one, beta, em, parry_orbit, parry_bits, nbits);

#define NBINS 403
	double mepa[NBINS];
	double mext[NBINS];
	double mexu[NBINS];

#define SCALE (4.0/3.0)
	for (int i=0; i< NBINS; i++)
	{
		double y = (((double) i) + 0.5) / ((double) NBINS);
		y *= SCALE;
		double pacc = 0.0;
		double xacc = 0.0;
		double uacc = 0.0;
		double bpn = 1.0;

		for (int k=0; k<20; k++)
		{
			double x = mpf_get_d(parry_orbit[k].get_mpf_t());
			if (y < x) pacc += bpn;

			for (int j=0; j<ntracks; j++)
			{
				std::vector<bool> bits = bitset[j];
				std::vector<mpf_class> orbit = orbit_set[j];
				std::vector<int> branch_points = branch_set[j];

				size_t nb = branch_points.size();
				size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
				if (orbit.size() <= norb) norb = orbit.size() -1;
				if ((int) norb < k) continue;

				std::vector<bool> ebits = ebitset[j];
				std::vector<mpf_class> erbit = erbit_set[j];

				double x = mpf_get_d(orbit[k].get_mpf_t());
				//double xu = mpf_get_d(erbit[k].get_mpf_t());

				if (bits[k] and y < x) xacc += bpn;
				//if (ebits[k] and y < xu) uacc += bpn;

#define ALMOST_WORKS_BUT_DOESNT
#ifdef ALMOST_WORKS_BUT_DOESNT

				// double scale = 1.0;
				bool branch = false;
				// int len = 0;
				for (int b = 0; b<(int)nb; b++)
					if (k == branch_points[b]) { branch = true; /*len=b;*/ break; }
				// if (branch) scale = -0.5/dbeta;

				// if (y < x) acc += bpn;
				// if (scale * y < scale * x) acc += bpn;
				if (not branch and y < x) uacc += bpn;
				if (branch)
				{
					if (bits[k] and y < x) uacc += bpn;
					if (not bits[k])
					{
						// if (y < x) uacc += bpn/dbeta;
						if (y < x) uacc += bpn;
						if (dbeta*y < x) uacc -= dbeta*bpn;
					}
				}
#endif // ALMOST_WORKS_BUT_DOESNT

			}

			bpn /= dbeta;
		}
		mepa[i] = pacc;
		mext[i] = xacc;
		mexu[i] = uacc;
	}

	double porm = 0.0;
	double xorm = 0.0;
	double uorm = 0.0;
	for (int i=0; i< NBINS; i++)
	{
		porm += mepa[i];
		xorm += mext[i];
		uorm += mexu[i];
	}

	for (int i=0; i< NBINS; i++)
	{
		mepa[i] *= NBINS/porm;
		mext[i] *= NBINS/xorm;
		mexu[i] *= NBINS/uorm;

		mepa[i] /= SCALE;
		mext[i] /= SCALE;
		mexu[i] /= SCALE;
	}

	for (int i=0; i< NBINS; i++)
	{
		double y = (((double) i) + 0.5) / ((double) NBINS);
		printf("%d	%g	%g	%g	%g\n", i, y*SCALE, mepa[i], mext[i], mexu[i]);
	}
}

// ================================================================
