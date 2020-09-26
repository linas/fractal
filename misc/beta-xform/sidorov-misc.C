/*
 * sidorov-misc.C
 * Bignum version of sidorov.C
 *
 * Understand the "gaps" aka the alternative beta-expansions
 * of any given number x, as described by N Sidorov.
 *
 * This includes graphing the gaps directly, and misc stuff
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#define HISTOGRAM_ORBITS
#include "sidorov-big.C"


#ifdef HISTOGRAM_ORBITS
int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K nbits\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	int nbits = atoi(argv[2]);

	do_init(nbits);

	mpf_class beta;
	make_random_bitsequence(beta, 2.0*Kay, nbits, NBINS);

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d nbits=%d\n#\n", Kay, em, nbits);

	// Where are the extended orbits going?
	// Draw a histogram
	double histo[NBINS];
	double tracknum[NBINS];
	double tracklen[NBINS];
	double tracklensq[NBINS];
	double tracklenrms[NBINS];
	double trackmax[NBINS];
	double xmean[NBINS];
	for (int i=0; i<NBINS; i++)
	{
		histo[i] = 0.0;
		histbase[i] = 0.0;

		xmean[i] = 0.0;
		tracknum[i] = 0.0;
		tracklen[i] = 0.0;
		tracklensq[i] = 0.0;
		tracklenrms[i] = 0.0;
		trackmax[i] = 0.0;
	}

#define MAXDEPTH 7
	// Distances to the branchpoints.
	double depthnum[MAXDEPTH];
	double depthsum[MAXDEPTH];
	for (int i=0; i<MAXDEPTH; i++)
	{
		depthnum[i] = 0.0;
		depthsum[i] = 0.0;
	}

	size_t tot_tracks = 0;
	size_t tot_tracklen = 0;
	size_t tot_tracklensq = 0;
	size_t tot_tracks_longest = 0;
	mpf_class ex;

#define NSAMP 512
	printf("# Sampled unit interval %d times\n#\n", NSAMP);

	for (int ibin=0; ibin<NBINS; ibin++)
	{
		// if (ibin%100 ==0) fprintf(stderr, "# orbits done %d of %d\n", ibin, NBINS);
		// fprintf(stderr, "# orbits done %d of %d\n", ibin, NBINS);
		double x = (((double) ibin) + 0.5)/ ((double) NBINS);

		double avg_x = 0.0;
		for (int nsamp=0; nsamp<NSAMP; nsamp++)
		{
			make_random_bitsequence(ex, x, nbits, NBINS);
			avg_x += mpf_get_d(ex.get_mpf_t());

			std::vector<std::vector<mpf_class>> orbit_set;
			std::vector<std::vector<bool>> bitset;
			std::vector<std::vector<int>> branch_set;
			beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

			tracknum[ibin] += bitset.size();
			tot_tracks += bitset.size();

			// Compute a histogram of the orbits. But do it
			// by summing only up to the last branch-point.
			// (else the greedy expansion will dominate).
			int longest = 0;
			int ntracks = bitset.size();
			for (int j=0; j<ntracks; j++)
			{
				std::vector<mpf_class> orbit = orbit_set[j];
				std::vector<int> branch_points = branch_set[j];
				int last = branch_points.back();
				tot_tracklen += last;
				tot_tracklensq += last * last;
				if (longest < last) longest = last;

				tracklen[ibin] += last;
				tracklensq[ibin] += last * last;
				if (trackmax[ibin] < longest) trackmax[ibin] = longest;

#define SCALE 1.3
				size_t nb = branch_points.size();
				size_t norb = 2*branch_points[nb-1] - branch_points[nb-2];
				if (orbit.size() <= norb) norb = orbit.size() -1;
				for (size_t k=1; k<=norb; k++)
				{
					double x = mpf_get_d(orbit[k].get_mpf_t());
					int bin = x * NBINS / SCALE;
					if (NBINS <= bin) bin=NBINS-1;
					histo[bin] += 1.0;
				}

				int nbranches = branch_points.size();
				for (int k=0; k<nbranches; k++)
				{
					depthnum[k] += 1.0;
					depthsum[k] += branch_points[k];
				}
			}
			tot_tracks_longest += longest;
		}

		avg_x /= NSAMP;
		xmean[ibin] = avg_x;

		tracklen[ibin] /= tracklen[ibin];
		tracklensq[ibin] /= tracklen[ibin];
		tracklenrms[ibin] = sqrt(tracklensq[ibin] - tracklen[ibin] * tracklen[ibin]);

	}
	double avg_tracks = ((double) tot_tracks) / (NBINS * NSAMP);
	double avg_tracklen = ((double) tot_tracklen) / tot_tracks;
	avg_tracklen *= log(2.0) / log(avg_tracks);

	double ms_tracklen = ((double) tot_tracklensq) / tot_tracks;
	ms_tracklen *= log(2.0) / log(avg_tracks);

	double rms_tracklen = sqrt(ms_tracklen - avg_tracklen*avg_tracklen);
	double avg_longest = tot_tracks_longest / (NBINS * NSAMP);

	printf("# Avg tracks/orbit: %g expect 2^%d=%d miss=%g avg tracklen: %g\n#\n",
	       avg_tracks, MAXDEPTH, 1<<MAXDEPTH, (1<<MAXDEPTH) - avg_tracks, avg_tracklen);
	printf("# RMS tracklen: %g Longest-avg: %g\n", rms_tracklen, avg_longest)'
	fprintf(stderr, "# Avg tracks/orbit: %g expect 2^%d=%d avg tracklen: %g\n",
	       avg_tracks, MAXDEPTH, 1<<MAXDEPTH, avg_tracklen);

// #define PRINT_HISTORGRAM
#ifdef PRINT_HISTORGRAM

	// Obtain a comparable number of counts for the baseline
	for (int i=0; i<NBINS; i++)
	{
		if (i%100 ==0) fprintf(stderr, "# baseline done %d of %d\n", i, NBINS);
		double x = (((double) i) + 0.5)/ ((double) NBINS);

#define BASE_SAMP 13000
// #define BASE_SAMP 100
		// 13000 seems to give a nice result...
		for (int j=0; j<BASE_SAMP; j++)
		{
			make_random_bitsequence(ex, x, nbits, NBINS);
			beta_sequence(ex, beta, em, nbits);
		}
	}

	// Normalize
	double cnt = 0.0;
	double bnt = 0.0;
	for (int i=0; i<NBINS; i++)
	{
		cnt += histo[i];
		bnt += histbase[i];
	}
	fprintf(stderr, "# histogram counts=%g baseline=%g\n", cnt, bnt);
	for (int i=0; i<NBINS; i++)
	{
		histo[i] *= NBINS/cnt;
		histbase[i] *= NBINS/bnt;
	}

	// Print the histogram
	double sum = 0.0;
	double bsu = 0.0;
	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		sum += histo[i] / NBINS;
		bsu += histbase[i] / NBINS;
		printf("%d	%g	%g %g	%g	%g %g	%g\n",
		       i, x*SCALE, histo[i], histbase[i], x, sum, bsu, tracklen[i]);
	}
#endif // PRINT_HISTOGRAM

#define PRINT_LENGTH_DISTRIBUTION
#ifdef PRINT_LENGTH_DISTRIBUTION
	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		printf("%d	%g	%g	%g	%g	%g\n",
		       i, x, xmean[i], tracklen[i], trackrms[i], trackmax[i]);
	}
#endif // PRINT_LENGTH_DISTRIBUTION

#ifdef PRINT_LENGTH
	// Normalize the distance to each branch point
	double prev = 0.0;
	for (int i=0; i<MAXDEPTH; i++)
	{
		depthsum[i] /= depthnum[i];
		depthnum[i] /= NBINS;

		printf("%d	%g	%g	%g\n", i, depthnum[i], depthsum[i], depthsum[i]-prev);

		prev = depthsum[i];
	}
#endif // PRINT_LENGTH
}
#endif // HISTOGRAM_ORBITS

// ================================================================
