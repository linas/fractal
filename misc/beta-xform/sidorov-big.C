/*
 * sidorov-big.C
 * Bignum version of sidorov.C
 *
 * Understand the "gaps" aka the alternative beta-expansions
 * of any given number x, as described by N Sidorov.
 *
 * This includes graphing the gaps directly, as well as graphing
 * the extended measure.
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#include <vector>

#include <gmpxx.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

gmp_randstate_t rstate;
mpf_class half;
mpf_class one;
mpf_class two;

void do_init(int nbits)
{
	gmp_randinit_default(rstate);
	mpf_set_default_prec(nbits);

	one = 1;
	two = 2;
	half = one / two;
}

// Given a double-precision value x, this will create a random
// bit-sequence that is nbits long, with the top 50 bits being
// those taken from the double-precision value x, and the rest
// randomly generated.  This is meant to provide a uniform sampling
// on the unit interval; equivalently, uniform sampling on the
// product space.
void make_random_bitsequence(mpf_class& val, double x, int nbits, int nbins)
{
	mpf_class tail;

	mpf_urandomb(tail.get_mpf_t(), rstate, nbits);

	// Rescale so that it's uniformly distributed over a bin.
	tail /= nbins;
	val = x;
	val += tail;
}

#define HISTOGRAM_ORBITS
#ifdef HISTOGRAM_ORBITS
	#define NBINS 403
	double histbase[NBINS];
#endif

// Compute the m from the sidorov paper. This is the length of the
// run of zeros we need to see, before exploring an alternate branch.
// It depends only on beta. Note that the Sidorov paper incorrectly
// states that the run of zeros is one less than m. This is false.
// The required run  of zeros has to be at least m.
// Note that:
// m=1 for K < 0.25 (1+sqrt(5)) = 0.80902
// m=2 up until about 0.878
// m=3 up until about 0.967
int emrun(double K)
{
	double beta = 2.0*K;
	double gold = 0.5 * (1.0 + sqrt(5));
	if (beta <= gold) return 1;

	double loga = (beta - 1.0) / (2.0-beta);
	loga = log(loga) / log(beta);
	loga = floor(loga) + 1.0;
	return (int) loga;
}

// ================================================================

// Generate the beta expansion in a greedy fashion.
// `y` is the number to expand.
// `K` is beta/2
// `nstart` is where to start writing the bit pattern
// `bits` is where to write the resulting bit-vector (it is never read)
// `orbit` is where to write the orbit sequence. (it is never read)
void greedy_expand(mpf_class y, mpf_class beta, int nstart, int nbits,
                   std::vector<mpf_class>& orbit,
                   std::vector<bool>& bits)
{
	// Generate the beta expansion bits in a greedy fashion.
	for (int i=nstart; i<nbits; i++)
	{
		orbit[i] = y;
		if (half <= y)
		{
			y -= half;
			bits[i] = 1;
		}
		else bits[i] = 0;
		y *= beta;
	}
}

// Recursively generate a collection of equivalent beta expansions at K.
// `y` is the number to expand.
// `K` is beta/2
// `em` is the number of zero bits to look for.
// `start` is where to start looking
// `gap` is the set to which the beta exapnsion is appended.
void beta_expand_rec(mpf_class y, mpf_class beta, int em, int start, int nbits,
                     std::vector<mpf_class> orbit,
                     std::vector<bool> greedy,
                     std::vector<int> branch_points,
                     int depth,
                     std::vector<std::vector<mpf_class>>& orbit_set,
                     std::vector<std::vector<bool>>& gap,
                     std::vector<std::vector<int>>& branch_set)
{
#define MAXDEPTH 10
	if (MAXDEPTH <= depth)
	{
		orbit_set.push_back(orbit);
		gap.push_back(greedy);
		branch_set.push_back(branch_points);
		return;
	}

	// Search for runs of length em.
	// The Sidorov paper has an error, the m is off by one.
	for (int i=start; i<nbits-em; i++)
	{
		if (1 == greedy[i])
		{
			bool found = true;
			for (int j=1; j<=em; j++)
			{
				if (1 == greedy[i+j])
					found = false;
			}
			if (found)
			{
				std::vector<bool> gapper = greedy;
				std::vector<mpf_class> lorbit = orbit;

				// Set to zero, and resume expansion.
				gapper[i] = 0;
				y = beta * orbit[i];
				lorbit[i] = y;

				branch_points.push_back(i);
				std::vector<int> lobran = branch_points;

				// Get the alternate expansion.
				greedy_expand(y, beta, i+1, nbits, lorbit, gapper);

				// Recurse
				beta_expand_rec(y, beta, em, i+1, nbits, orbit, greedy,
				                branch_points, depth+1, orbit_set, gap, branch_set);
				beta_expand_rec(y, beta, em, i+1, nbits, lorbit, gapper,
				                lobran, depth+1, orbit_set, gap, branch_set);
				return;
			}
		}
	}
}

void beta_expand(mpf_class y, mpf_class beta, int em,
                 std::vector<std::vector<mpf_class>>& orbit_set,
                 std::vector<std::vector<bool>>& gap,
                 std::vector<std::vector<int>>& branch_set,
                 int nbits)
{
	std::vector<bool> bitseq;
	std::vector<mpf_class> orbit;
	std::vector<int> branch_points;
	bitseq.resize(nbits);
	orbit.resize(nbits);

	for (int i=0; i< nbits; i++)
	{
		orbit[i] = 0;
	}

	// First, get the baseline (greedy) orbit.
	greedy_expand(y, beta, 0, nbits, orbit, bitseq);

	// And now recursively get the rest.
	beta_expand_rec(y, beta, em, 0, nbits, orbit, bitseq, branch_points,
	                0 /* depth*/, orbit_set, gap, branch_set);

#if 0 // DEBUG
	// Debug print
	for (size_t n=0; n<gap.size(); n++)
	{
		std::vector<bool> bits = gap[n];
		printf("# n=%lu ", n);
		for (int i=0; i<NBITS; i++)
			printf("%d", (int) bits[i]);
		printf("\n");
	}
	printf("#\n#\n");
#endif
}

std::vector<bool> beta_sequence(mpf_class y, mpf_class beta,
                                int em, int nbits)
{
	std::vector<bool> bitseq;
	std::vector<mpf_class> orbit;
	bitseq.resize(nbits);
	orbit.resize(nbits);

	for (int i=0; i< nbits; i++)
	{
		orbit[i] = 0;
	}

	// Get the baseline (greedy) orbit.
	greedy_expand(y, beta, 0, nbits, orbit, bitseq);

#ifdef HISTOGRAM_ORBITS
#define SCALE 0.75
	for (int i=0; i< nbits; i++)
	{
		double x = mpf_get_d(orbit[i].get_mpf_t());
		int bin = x * NBINS * SCALE;
		if (NBINS <= bin) bin=NBINS-1;
		histbase[bin] += 1.0;
	}
#endif // HISTOGRAM_ORBITS

	return bitseq;
}

// ================================================================

#if 0
// Sum the bit-sequence, returning the sum.
// This returns 0.5 * sum_i=0 b[i] (2J)^-i
//
double beta_sum(std::vector<bool> bits, double Jay)
{
	double acc = 1.0e-30;
	for (int i=0; i<NBITS; i++)
	{
		acc *= 1.0 / (2.0*Jay);
		if (bits[NBITS-i-1])
		{
			acc += 0.5;
		}
	}
	return acc;
}
#endif

// ================================================================

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

#ifdef HISTOGRAM_ORBITS
	// Where are the extended orbits going?
	// Draw a histogram
	double histo[NBINS];
	double tracklen[NBINS];
	for (int i=0; i<NBINS; i++)
	{
		histo[i] = 0.0;
		histbase[i] = 0.0;
		tracklen[i] = 0.0;
	}

	// Distances to the branchpoints.
	double tracknum[MAXDEPTH+2];
	double tracksum[MAXDEPTH+2];
	for (int i=0; i<MAXDEPTH+2; i++)
	{
		tracknum[i] = 0.0;
		tracksum[i] = 0.0;
	}

	size_t tot_tracks = 0;
	size_t tot_tracklen = 0;
	mpf_class ex;
	for (int ibin=0; ibin<NBINS; ibin++)
	{
		if (ibin%100 ==0) fprintf(stderr, "# orbits done %d of %d\n", ibin, NBINS);
		double x = (((double) ibin) + 0.5)/ ((double) NBINS);
		make_random_bitsequence(ex, x, nbits, NBINS);

		std::vector<std::vector<mpf_class>> orbit_set;
		std::vector<std::vector<bool>> bitset;
		std::vector<std::vector<int>> branch_set;
		beta_expand(ex, beta, em, orbit_set, bitset, branch_set, nbits);

		tot_tracks += bitset.size();

		// Compute a histogram of the orbits. But do it only by 
		// summing up to the last branch-point.
		int ntracks = bitset.size();
		for (int j=0; j<ntracks; j++)
		{
			std::vector<mpf_class> orbit = orbit_set[j];
			std::vector<int> branch_points = branch_set[j];
			int last = branch_points.back();
			tot_tracklen += last;
			tracklen[ibin] += last;
			for (int k=0; k<=last+1; k++)
			{
				double x = mpf_get_d(orbit[k].get_mpf_t());
				int bin = x * NBINS * SCALE;
				if (NBINS <= bin) bin=NBINS-1;
				histo[bin] += 1.0;
			}

			int nbranches = branch_points.size();
			for (int k=0; k<nbranches; k++)
			{
				tracknum[k] += 1.0;
				tracksum[k] += branch_points[k];
			}
		}
		tracklen[ibin] /= ntracks;
	}
	double avg_tracks = ((double) tot_tracks) / NBINS;
	double avg_tracklen = ((double) tot_tracklen) / tot_tracks;
	printf("# Avg tracks/orbit: %g avg tracklen: %g\n",
	       avg_tracks, avg_tracklen);
	fprintf(stderr, "# Avg tracks/orbit: %g avg tracklen: %g\n",
	       avg_tracks, avg_tracklen);

#ifdef PRINT_HISTORGRAM
	int navg = avg_tracks;

	// Obtain a comparable number of counts for the baseline
	for (int i=0; i<NBINS; i++)
	{
		if (i%100 ==0) fprintf(stderr, "# baseline done %d of %d\n", i, NBINS);
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		for (int j=0; j<navg; j++)
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
	fprintf(stderr, "# histogram counts=%d baseline=%d\n", (int) cnt, (int) bnt);
	for (int i=0; i<NBINS; i++)
	{
		histo[i] *= NBINS/cnt;
		histbase[i] *= NBINS/bnt;
	}

	// Print the histogram
	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		x /= SCALE;
		printf("%d	%g	%g	%g	%g\n", i, x, histo[i], histbase[i], tracklen[i]);
	}
#endif // PRINT_HISTOGRAM

	// Normalize the distance to each branch point
	for (int i=0; i<MAXDEPTH+2; i++)
	{
		tracksum[i] /= tracknum[i];
		tracknum[i] /= NBINS;
	}

	// And now the sum from the paper
	double pee[MAXDEPTH+2];
	double peel = 0.0;
	for (int i=0; i<MAXDEPTH+2; i++)
	{
		peel += tracksum[i] * tracknum[i];
		pee[i] = peel / tracknum[i];

		printf("%d	%g	%g	%g\n", i, tracknum[i], tracksum[i], pee[i]);
	}
#endif
}

// ================================================================
