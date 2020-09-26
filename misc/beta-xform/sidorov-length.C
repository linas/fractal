/*
 * sidorov-length.C
 * Compute distance betweeen branch points.
 *
 * Linas Vepstas Sept 2020
 */

#include "sidorov-big.C"


int main (int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, "Usage: %s nbits depth npts\n", argv[0]);
		exit (1);
	}
	int nbits = atoi(argv[1]);
	int maxdepth = atoi(argv[2]);
	int npts = atoi(argv[3]);

	do_init(nbits);

#define MAXDEPTH maxdepth
#define NSAMP 200
	printf("#\n# Dataset for average track length as function of K\n");
	printf("#\n# Computed for %d/(beta-1)(2-beta) bits precision\n", nbits);
	printf("# Sampled unit interval %d times\n", NSAMP);
	printf("#\n# Column labels:\n");
	printf("# kay, avg-tracks/orbit expect 2^%d=%d, deficit, avg-tracklen, avg-longest, rms-len\n#\n",
		MAXDEPTH, 1<<MAXDEPTH);

	for (int i=0; i<npts; i++)
	{
		double dbeta = 1.0 + (((double) i) + 0.5)/ ((double) npts);

		mpf_class beta;
		make_random_bitsequence(beta, dbeta, nbits, npts);
		int em = emrun(dbeta/2.0);

		fprintf(stderr, "working beta=%g\n",  dbeta);

		double tot_tracks = 0.0;
		double tot_tracklen = 0.0;
		double tot_tracklensq = 0.0;
		double tot_tracks_longest = 0.0;

		double avg_x = 0.0;
		for (int nsamp=0; nsamp<NSAMP; nsamp++)
		{
			mpf_class ex;
			make_random_bitsequence(ex, 0.0, nbits, 1);
			avg_x += mpf_get_d(ex.get_mpf_t());

			std::vector<std::vector<mpf_class>> orbit_set;
			std::vector<std::vector<bool>> bitset;
			std::vector<std::vector<int>> branch_set;
			beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

			tot_tracks += bitset.size();

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
			}
			tot_tracks_longest += longest;
		}

		avg_x /= NSAMP;
printf("duuude avgx=%g\n", avg_x);

		// Collect up tracklen stats.
		double avg_tracks = ((double) tot_tracks) / NSAMP;
		double avg_tracklen = ((double) tot_tracklen) / tot_tracks;
		avg_tracklen *= log(2.0) / log(avg_tracks); 

		double ms_tracklen = ((double) tot_tracklensq) / tot_tracks;
		ms_tracklen *= log(2.0) / log(avg_tracks); 

		double rms_tracklen = sqrt(ms_tracklen - avg_tracklen*avg_tracklen);

		double avg_longest = tot_tracks_longest / NSAMP;

		printf("%g	%g	%g	%g	%g	%g\n", 0.5*dbeta,
		       avg_tracks, (1<<MAXDEPTH) - avg_tracks, avg_tracklen, avg_longest, rms_tracklen);
		fflush(stdout);
	}
}


// ================================================================
