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

	int em = emrun(Kay);
	double dbeta = 2.0*Kay;
	double dalpha = 0.5 * (1.0 + pow(dbeta, -em));

	printf("#\n# K=%g beta=%g alpha=%g m=%d nbits=%d\n#\n", Kay, dbeta, dalpha, em, nbits);

#define MAXDEPTH maxdepth

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
	printf("# ntracks=%d etrack=%d\n", ntracks, etracks);

	std::vector<mpf_class> parry_orbit;
	std::vector<bool> parry_bits;
	beta_sequence(one, beta, em, parry_orbit, parry_bits, nbits);

	// Integral of the Parry measure.
	double pgral = 0.0;
	double bpn = 1.0;
	for (int k=0; k< (int) parry_orbit.size(); k++)
	{
		double x = mpf_get_d(parry_orbit[k].get_mpf_t());
		pgral += x*bpn;
// printf("# %d %g %g %g\n", k, x, bpn, pgral);
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
		double apn = 1.0;
		double cpn = 1.0;

		for (int k=0; k<20; k++)
		{
			double t = mpf_get_d(parry_orbit[k].get_mpf_t());
			if (y < t) pacc += bpn;

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
				if (y < xu) uacc += bpn;    // works best for beta> 1.5
				// if (y < xu/dbeta) uacc += apn;  // works best for beta<1.5
				// if (y < xu and not bits[k]) uacc += bpn; // works best for beta<1.5

				// fails
				// if (y < xu) dacc += 2*(bits[k]-0.5) * bpn; // fail

				if (y < x) dacc += bpn;
				// if (y < xu) dacc += bpn;
				// if (y < xu) dacc += (dbeta-1.0)*bpn;
				// if (y < xu) dacc += (1.0-dalpha)*bpn;

				// if (y < xu - 0.5*dbeta) dacc -= bpn;

				// Wow, this one is almost perfect!
				// if (y < xu - x) dacc -= bpn;

				// Wow! This is even better! ... but still off
				// if (y < xu - x) dacc -= bpn /  (2.0*dalpha);

				// if (y < xu - x) dacc -= bpn * (0.5 * dbeta) / (2.0*dalpha);

				if (y < xu - x) dacc -= bpn * 0.25* (0.5 * dbeta) / (2.0*dalpha);

#ifdef CLOSE_AT_BETA_12
				// The pure if (y < x) dacc += bpn; fits great at x<0.25
				// The pure if (y < xu) dacc += bpn; fits great at x=0.5
				// The if (y < xu - x) is needed for the bump.
				// the assembly of them all .. meh.
				if (y < x) dacc += bpn;
				if (y < xu - x) dacc -= bpn * 0.25* (0.5 * dbeta) / (2.0*dalpha);
#endif // CLOSE_AT_BETA_12

#ifdef MORE_OR_LESS_PERFECT_AT_BETA_188
				// This seems to give a more-or-less perfect fit
				// at beta=1.88 with maybe imperfection above x>0.94 ??
				// Argh ... but fails once moving away from there...
				if (y < x) dacc += bpn;
				if (y < xu) dacc += (dbeta-1.0)*bpn;
				if (y < xu - x) dacc -= bpn * (0.5 * dbeta) / (2.0*dalpha);
#endif // MORE_OR_LESS_PERFECT_AT_BETA_1.88



// #define ALMOST_WORKS_BUT_DOESNT
#ifdef ALMOST_WORKS_BUT_DOESNT

// Should be titled: "almost works but can't".  The clue is given
// at Kay=0.94 (beta=1.88) where the first branch occurs at k=47
// and there are no earlier branches, and so the code below is effectively
// never entered. Yet, clearly, the density has whack steps. So the
// fact that the below "almost" worked really is accidental.

// Note that x < xu until after the first branch point
// so that is kind-of another way of detecting branch points.
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

						// These almost work because they shift big step from
						// 0.75 to 0.5 (when beta=1.5) but then there's to much fiddle
						// and for beta=0.85 it just steps wrong...
						//if (dbeta*y < x) dacc -= bpn ;// dbeta; // almost; but that step
						// if (y < x/dbeta) dacc -= bpn;

						// This is almost perfect for beta=1.5
						// Like above, the right place for the step, and without grunge!
						// if (y < 0.5) dacc -= bpn /dbeta ;

						// if (y < xu) dacc += bpn ;//dbeta;  // sometimes OK adjust

						// strangely, this is wrong under the branch...
						// which is weird given later ones work for beta=1.7
						// if (y < x/dbeta) dacc -= bpn ; //*dbeta;
						// ---------------------


						// This alone gives a nearly perfect fit for for beta=1.5
						// for 0 < y < 0.75 and comes close for the rest.
						// It's overkill for beta=1.3
						// if (y < 0.5) dacc -=  apn/(dbeta*dalpha);

						// if (y < 0.5) dacc -=  apn/(dbeta*dalpha);
						// xxx if (y < 0.5) dacc -= bpn/(2*dbeta);

						// This is almost identical to above, but not quite.
						// This is better for small beta, but still overkill
						// if (y < 0.5) dacc -= dbeta*cpn;

						// These two combined give a nearly perfect fit for beta=1.7
						// for y greater than 0.5 and come close for the rest.
						// if (y < 0.5) dacc -=  apn/(dbeta*dalpha);
						// if (y < 0.5/dbeta) dacc -= apn/dbeta ;

						// These three combined give an OK fit for beta=1.76 to 1.86
						// So there is an m=3 effect at work, here.
						// if (y < 0.5) dacc -=  apn/(dbeta*dalpha);
						// if (y < 0.5/dbeta) dacc -= apn/dbeta ;
						// if (y < 0.5/(dbeta*dbeta)) dacc -= apn ;

						// These four combined give an OK fit for beta=1.87 to 1.93
						// So there is an m=4 effect at work, here.
						// XXX FAIL!
						if (y < 0.5) dacc -=  apn/(dbeta*dalpha);
						if (y < 0.5/dbeta) dacc -= apn/dbeta ;
						if (y < 0.5/(dbeta*dbeta)) dacc -= apn ;
						if (y < 0.5/(dbeta*dbeta*dbeta)) dacc -= apn ; // ???

						// Meh.
						//if (y < 0.5) dacc -= dbeta*cpn;
						//if (y < 0.5/dbeta) dacc -= dalpha*dbeta*cpn ;
					}

					// this fails for beta=1.5 but needed for beta=1.7
					// also fine details are wrong; below is better.
					// if (y < 0.5*x) dacc -= bpn *dbeta;

					// this fails for beta=1.5 but needed for beta=1.7
					// This suggests that its actual length that matters.
					// if (y < 0.5/dbeta) dacc -= bpn * dbeta;
					// if (y < 0.5/dbeta) dacc -= bpn /dbeta;


					// if (y < 0.5) dacc -=  bpn /(dbeta*dalpha) ;
					// if (y < 0.5) dacc -=  apn/dbeta;
				}

				// Almost perfect for beta=1.5 when taken outside of the if's
				// The fact that three divides are required hints that
				// its actually the length at play, cause length=3.27 at beta=1.5
				// and since its the length, it must need the branch.
				// if (y < 0.5) dacc -= bpn /dbeta/dbeta/dbeta ;

				// This plus above is almost perfect for beta=1.7
				// Again this hints at length, due to length change here...
				//if (y < 0.5/dbeta) dacc -= bpn /dbeta/dbeta;
#endif // ALMOST_WORKS_BUT_DOESNT

			}

			bpn /= dbeta;
			apn /= dbeta*dalpha;
			cpn /= 2.0*dbeta*dalpha;
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

	double ow = NBINS/SCALE;
	for (int i=0; i< NBINS; i++)
	{
		mepa[i] *= ow/porm;
		meda[i] *= ow/dorm;
		mext[i] *= ow/xorm;
		mexu[i] *= ow/uorm;

mepa[i] /= dbeta;
	}

	printf("# Parry meas norm=%g\n", porm/ow);

	for (int i=0; i< NBINS; i++)
	{
		double y = (((double) i) + 0.5) / ((double) NBINS);
		printf("%d	%g	%g	%g	%g	%g\n", i, y*SCALE, mepa[i], meda[i], mext[i], mexu[i]);
	}
}

// ================================================================
