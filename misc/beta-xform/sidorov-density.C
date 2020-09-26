/*
 * sidorov-density.C
 * Graph the extened density.
 *
 * Linas Vepstas Sept 2020
 */

#include <pthread.h>
#include "brat.h"

#include "sidorov-big.C"


static void extended_measure (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double Kay,
                             int itermax,
                             double param)
{
	static mpf_class beta;
	static std::vector<double> Kvec;
	static std::vector<double> trackvec;
	int nbits = itermax;

// #define NBINS array_size
#define NBINS 67
	static bool init=false;
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	if (not init)
	{
		pthread_mutex_lock(&mutex);
		if (not init)
		{
			do_init(nbits);
// #define NSAMP 16
#define NSAMP 3
#define MAXDEPTH 16
			printf("#\n# Dataset for average track length as function of K\n");
			printf("#\n# Computed for %d bits precision\n", nbits);
			printf("# Sampled unit interval %d times %d bins\n", NSAMP, NBINS);
			printf("#\n# Column labels:\n");
			printf("# kay, avg-tracks/orbit expect 2^%d=%d, deficit, avg-tracklen, avg-longest, rms-len\n#\n",
				MAXDEPTH, 1<<MAXDEPTH);


			init = true;
		}
		pthread_mutex_unlock(&mutex);
	}

	make_random_bitsequence(beta, 2.0*Kay, nbits, NBINS);
	int em = emrun(Kay);

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	fprintf(stderr, "working K=%g\n",  Kay);

	double tot_tracks = 0.0;
	double tot_tracklen = 0.0;
	double tot_tracklensq = 0.0;
	double tot_tracks_longest = 0.0;
	double xmean[NBINS];

	mpf_class ex;
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
				tot_tracklensq += last*last;
				if (longest < last) longest = last;
#define SCALE 1.3
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
			tot_tracks_longest += longest;
		}

		avg_x /= NSAMP;
		xmean[ibin] = avg_x;
	}

	// Normalize
	double nobs = 0;
	for (int j=0; j<NBINS; j++)
		nobs += array[j];
	for (int j=0; j<NBINS; j++)
		array[j] *= NBINS/nobs;


	// Collect up tracklen stats.
	double avg_tracks = ((double) tot_tracks) / (NBINS * NSAMP);
	double avg_tracklen = ((double) tot_tracklen) / tot_tracks;
	avg_tracklen *= log(2.0) / log(avg_tracks); 

	double ms_tracklen = ((double) tot_tracklensq) / tot_tracks;
	ms_tracklen *= log(2.0) / log(avg_tracks); 

	double rms_tracklen = sqrt(ms_tracklen - avg_tracklen*avg_tracklen);

	double avg_longest = tot_tracks_longest / (NBINS * NSAMP);

	printf("%g	%g	%g	%g	%g	%g\n", Kay,
	       avg_tracks, (1<<MAXDEPTH) - avg_tracks, avg_tracklen, avg_longest, rms_tracklen);
	fflush(stdout);
}

DECL_MAKE_BIFUR(extended_measure)

// ================================================================
