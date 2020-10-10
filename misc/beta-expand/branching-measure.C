/*
 * branching-measure.C
 * Try to guess what the branching measure is.
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

	// Orbits of 1.0 aka orbits of "lower bound" beta/2
	// That is, we do the first iteration by hand.
	std::vector<std::vector<mpf_class>> orbit_set;
	std::vector<std::vector<bool>> bitset;
	std::vector<std::vector<int>> branch_set;
	std::vector<std::vector<bool>> gamma_set;
	mpf_class one = Kay;
	mpf_class beta = dbeta;
	beta_expand(one, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, gamma_set, nbits);

	// Orbits of the upper bound, alpha*beta = beta/2 (1+\beta^-m)
	std::vector<std::vector<mpf_class>> erbit_set;
	std::vector<std::vector<bool>> ebitset;
	std::vector<std::vector<int>> ebr_set;
	std::vector<std::vector<bool>> egam_set;

	mpf_class alpha = Kay * (1.0 + pow(dbeta, -em));
	beta_expand(alpha, beta, em, MAXDEPTH, erbit_set, ebitset, ebr_set, egam_set, nbits);

	int ntracks = bitset.size();
	int etracks = ebitset.size();
	printf("# ntracks=%d etrack=%d\n", ntracks, etracks);

	std::vector<mpf_class> parry_orbit;
	std::vector<bool> parry_bits;
	beta_sequence(one, beta, em, parry_orbit, parry_bits, nbits);

	// Integral of the Parry measure.
	double pgral = 0.0;
	double bpn = 0.5 * dbeta;
	for (int k=0; k< (int) parry_orbit.size(); k++)
	{
		double x = mpf_get_d(parry_orbit[k].get_mpf_t());
		pgral += x*bpn;
		bpn /= dbeta;
	}
	printf("# Parry measure integral=%g\n", pgral);

#define NBINS 1803
	double mepa[NBINS];
	double meda[NBINS];
	double mext[NBINS];
	double mexu[NBINS];

#define SCALE (4.0/3.0)
	for (int i=0; i< NBINS; i++)
	{
		double y = (((double) i) + 0.5) / ((double) NBINS);
		y *= SCALE;
		double pacc = 0.0;
		double dacc = 0.0;
		double xacc = 0.0;
		double uacc = 0.0;
		double bpn = 1.0;

		for (int k=0; k<20; k++)
		{
			double t = mpf_get_d(parry_orbit[k].get_mpf_t());
			if (y < t) pacc += bpn;

			for (int j=0; j<ntracks; j++)
			{
				std::vector<bool>& bits = bitset[j];
				std::vector<mpf_class>& orbit = orbit_set[j];
				std::vector<int>& branch_points = branch_set[j];

				size_t nb = branch_points.size();
				size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
				if (orbit.size() <= norb) norb = orbit.size() -1;
				if ((int) norb < k) continue;

				// std::vector<bool>& ebits = ebitset[j];
				std::vector<mpf_class>& erbit = erbit_set[j];
				std::vector<int>& ebranches = branch_set[j];

				double x = mpf_get_d(orbit[k].get_mpf_t());
				double xu = mpf_get_d(erbit[k].get_mpf_t());

				if (y < x) xacc += bpn;
				if (y < xu) uacc += bpn;

#define ALMOST_WORKS_BUT_DOESNT
#ifdef ALMOST_WORKS_BUT_DOESNT

				bool branch = false;
				int b;
				for (b = 0; b<(int)nb; b++)
					if (k == branch_points[b]) { branch = true; break; }

				// We expect branch points to e in the same locations
				// for upper and lower iterates, and that does seem to hold.
				if (branch and k != ebranches[b])
				{
					printf("branch fail! \n");
					exit(1);
				}

				if (y < x) dacc += bpn;
				if (y < xu) dacc += bpn;

				if (branch)
				{
					// if (y < x) dacc -= bpn;

					if (not bits[k])
					{
					// if (dbeta*y < x) dacc -= bpn / dbeta; // almost

					// if (x < y) dacc -= bpn / dbeta; // Nooo
					// if (x < y and y < dbeta*xu) dacc += bpn; // Noo

					// if (x<y and y < xu) dacc += bpn;
					// if (y < xu) dacc += bpn /dbeta;
					// if (y < x) dacc -= bpn;  // Nooo
					// if (y < xu) dacc -= bpn /dbeta; // too negative

					// if (dbeta*y < x) dacc -= bpn / dbeta; // almost; but that step
					// if (Kay*y < x) dacc -= bpn; // ok but No, whack.
					// if (y < xu) dacc += bpn /dbeta;

					//if (dbeta*y < x) dacc -= bpn ;// dbeta; // almost; but that step
					//if (y < xu) dacc += bpn ;//dbeta;
					if (y < 0.5) dacc -= bpn *dbeta;
					}

				}
#endif // ALMOST_WORKS_BUT_DOESNT

			}

			bpn /= dbeta;
		}
		mepa[i] = pacc;
		meda[i] = dacc;
		mext[i] = xacc;
		mexu[i] = uacc;
	}

	double porm = 0.0;
	double dorm = 0.0;
	double xorm = 0.0;
	double uorm = 0.0;
	for (int i=0; i< NBINS; i++)
	{
		porm += mepa[i];
		dorm += meda[i];
		xorm += mext[i];
		uorm += mexu[i];
	}

	for (int i=0; i< NBINS; i++)
	{
		mepa[i] *= NBINS/porm;
		meda[i] *= NBINS/dorm;
		mext[i] *= NBINS/xorm;
		mexu[i] *= NBINS/uorm;
	}

	printf("# Parry meas norm=%g\n", porm/NBINS);

	for (int i=0; i< NBINS; i++)
	{
		double y = (((double) i) + 0.5) / ((double) NBINS);
		printf("%d	%g	%g	%g	%g	%g\n", i, y*SCALE, mepa[i], meda[i], mext[i], mexu[i]);
	}
}

// ================================================================
