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

#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

gmp_randstate_t rstate;
mpf_t half;
mpf_t one;
mpf_t two;

void do_init(int nbits)
{
	gmp_randinit_default(rstate);
	mpf_set_default_prec(nbits);

	mpf_init(one);
	mpf_init(half);
	mpf_init(two);
	mpf_set_ui(one, 1);
	mpf_set_ui(two, 2);
	mpf_div_ui(half, one, 2);
}

// Given a double-precision value x, this will create a random
// bit-sequence that is nbits long, with the top 50 bits being
// those taken from the double-precision value x, and the rest
// randomly generated.  This is meant to provide a uniform sampling
// on the unit interval; equivalently, uniform sampling on the
// product space.
void make_random_bitsequence(mpf_t& val, double x, int nbits)
{
	mpf_t tail;
	mpf_init(tail);

	mpf_set_d(val, x);
	mpf_urandomb(tail, rstate, nbits);

	// Keep the top 12 decimal digits of x
	unsigned long digs = 1000000;
	digs *= 1000000;
	mpf_div_ui(tail, tail, digs);
	mpf_add(val, val, tail);
}

#define NBINS 403
double histo[NBINS];
double histbase[NBINS];

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
void greedy_expand(double y, double K, int nstart,
                   std::vector<double>& orbit,
                   std::vector<bool>& bits)
{
	// Generate the beta expansion bits in a greedy fashion.
	for (int i=nstart; i<NBITS; i++)
	{
		orbit[i] = y;
		if (0.5 <= y)
		{
			y -= 0.5;
			bits[i] = 1;
		}
		else bits[i] = 0;
		y *= 2.0*K;
	}
}

// Recursively generate a collection of equivalent beta expansions at K.
// `y` is the number to expand.
// `K` is beta/2
// `em` is the number of zero bits to look for.
// `start` is where to start looking
// `gap` is the set to which the beta exapnsion is appended.
void beta_expand_rec(double y, double K, int em, int start,
                     std::vector<double> orbit,
                     std::vector<bool> greedy,
                     int depth,
                     std::vector<std::vector<bool>>& gap)
{
#define MAXDEPTH 8
	if (MAXDEPTH < depth) return;

#define HISTOGRAM_ORBITS
#ifdef HISTOGRAM_ORBITS
#define SCALE 0.75
	for (int i=0; i< NBITS; i++)
	{
		double x = orbit[i];
		int bin = x * NBINS * SCALE;
		if (NBINS <= bin) bin=NBINS-1;
		histo[bin] += 1.0;
	}
#endif

	// Search for runs of length em.
	// The Sidorov paper has an error, the m is off by one.
	for (int i=start; i<NBITS-em; i++)
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
				std::vector<double> lorbit = orbit;

				// Set to zero, and resume expansion.
				gapper[i] = 0;
				y = 2.0*K * orbit[i];
				lorbit[i] = y;

				// Get the alternate expansion.
				greedy_expand(y, K, i+1, lorbit, gapper);
				gap.push_back(gapper);

				// Recurse
				beta_expand_rec(y, K, em, i+1, orbit, greedy, depth+1, gap);
				beta_expand_rec(y, K, em, i+1, lorbit, gapper, depth+1, gap);
				return;
			}
		}
	}
}

std::vector<std::vector<bool>> beta_expand(double y, double K, int em)
{
	std::vector<std::vector<bool>> gap;
	std::vector<bool> bits;
	std::vector<double> orbit;
	bits.resize(NBITS);
	orbit.resize(NBITS);

	// First, get the baseline (greedy) orbit.
	greedy_expand(y, K, 0, orbit, bits);
	gap.push_back(bits);

#ifdef HISTOGRAM_ORBITS
	for (int i=0; i< NBITS; i++)
	{
		double x = orbit[i];
		int bin = x * NBINS * SCALE;
		if (NBINS <= bin) bin=NBINS-1;
		histbase[bin] += 1.0;
	}
#endif // HISTOGRAM_ORBITS

	// And now recursively get the rest.
	beta_expand_rec(y, K, em, 0, orbit, bits, 0, gap);

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
	return gap;
}

// ================================================================

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

// ================================================================

#if 1 // BRINGUP_DEBUG
// Basic unit-test - debugging code
int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K x\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double x = atof(argv[2]);

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d x=%g\n#\n", Kay, em, x);

#if CHECK_GROUND_TRUTH
	// Perform the basic sanity check that everything is OK.
	int npts = 313;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double z = sdr (x, Kay, em);
		printf("%d	%g	%g\n", i, x, z);
	}
#endif

#ifdef CHECK_MORE
	// More code validation
	int npts = 313;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);

		std::vector<std::vector<bool>> bitset;
		bitset = beta_expand(x, Kay, em);

		std::vector<bool> bits = bitset[0];
		double y = beta_sum(bits, Kay);
		printf("%d	%g	%g", i, x, y);
		for (size_t n=1; n<bitset.size(); n++)
		{
			std::vector<bool> bits = bitset[n];
			double z = beta_sum(bits, Kay);
			if (1.0e-7 < fabs(y-z))
				printf(" FAIL: %g", y-z);
		}
		printf("\n");
	}
#endif

#if GRAPH_GAPS
	std::vector<std::vector<bool>> bitset;
	bitset = beta_expand(x, Kay, em);

	int npts = 313;
	for (int i=0; i<npts; i++)
	{
		double Jay = (((double) i) + 0.5)/ ((double) npts);
		Jay = Kay + (1.0-Kay)*Jay;

		printf("%d	%g", i, Jay);
		for (size_t n=0; n<bitset.size(); n++)
		{
			std::vector<bool> bits = bitset[n];
			double y = beta_sum(bits, Jay);
			printf("	%g", y);
		}
		printf("\n");
	}
#endif

#define EXTENDED_MEASURE
#ifdef EXTENDED_MEASURE
	// Where are the extended orbits going?
	// Draw a histogram
	for (int i=0; i<NBINS; i++)
	{
		histo[i] = 0.0;
		histbase[i] = 0.0;
	}

	for (int i=0; i<NBINS; i++)
	{
		// if (i%100 ==0) printf("# done %d of %d\n", i, NBINS);
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		beta_expand(x, Kay, em);
	}

	// Normalize
	double cnt = 0.0;
	double bnt = 0.0;
	for (int i=0; i<NBINS; i++)
	{
		cnt += histo[i];
		bnt += histbase[i];
	}
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
		printf("%d	%g	%g	%g\n", i, x, histo[i], histbase[i]);
	}
#endif
}
#endif

// ================================================================

#ifdef DRAW_THE_GAPS

// This block of the code used to draw the sidorov "gaps"
// in the paper. Full collor. Tne actual paper used.
//  ./sidorov sid 800 800 1 0.5 0.0 1 0.75
//

#include "brat.h"
void MakeHisto (char * name,
                float *array,
                int      sizex,
                int      sizey,
                double x_center,
                double y_center,
                double width,
                double height,
                int itermax,
                double Kay)
{
	// double beta = 2.0*Kay;
	int em = emrun(Kay);

	printf("Tongues for K=%g em=%d\n", Kay, em);
	int tot_tracks = 0;
	for (int i=0; i<sizex; i++)
	{
		if (0 == i%100) printf("Start working column %d\n", i);
		double x = ((double) i + 0.5) / ((double) sizex);

		std::vector<std::vector<bool>> bitset;
		bitset = beta_expand(x, Kay, em);
		double rtracks = 1.0 / ((double) bitset.size());
		tot_tracks += bitset.size();

		for (int j=0; j<sizey; j++)
		{
			double y = ((double) j + 0.5) / ((double) sizey);
			double Jay = Kay + (1.0 - Kay) * y;
			for (size_t n=0; n<bitset.size(); n++)
			{
				std::vector<bool> bits = bitset[n];
				double p = beta_sum(bits, Jay);
				int ni = sizex * p;

				// frac runs from 0 to 1 from left of pixel to right.
				double frac = sizex * p - ni;

				// ne is the pixel to the left or to the right.
				int ne = ni;
				if (frac < 0.5) ne--; else ne++;

				if (sizex <= ni) ni = sizex - 1;
				if (sizex <= ne) ne = sizex - 1;
				if (ne < 0) ne = 0;

				// weight runs from 1.0 in the center of the pixel
				// to 0.5 at eigher edge of the pixel
				double weight = 1.0 - fabs(frac - 0.5);

				// if frac = 0.5 i.e. dead center, the center pixel
				// gets all the weight. Else the neighbor gets some.
				array[j*sizex + ni] += weight * rtracks;
				array[j*sizex + ne] += (1.0-weight) * rtracks;
			}
		}
	}
	double avg = ((double) tot_tracks) / ((double) sizex);
	printf("Graph shows an average of %g tracks per expansion\n", avg);
}

#endif
