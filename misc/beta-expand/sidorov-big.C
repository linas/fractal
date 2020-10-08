/*
 * sidorov-big.C
 * Bignum version of sidorov.C
 *
 * Understand the "gaps" aka the alternative beta-expansions
 * of any given number x, as described by N Sidorov.
 * Contains only the core recursive code.
 *
 * Linas Vepstas Sept 2020
 */

#include <vector>

#include <gmpxx.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "emrun.C"

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
// `greedy` is the bit expansion
// `orbit` is the collection of tee values
// `branch_points` locates where the branches are
// `gamm` is the left-right branching decisions
static
void beta_expand_rec(mpf_class y, mpf_class beta, int em, int start, int nbits,
                     std::vector<mpf_class> orbit,
                     std::vector<bool> greedy,
                     std::vector<int> branch_points,
                     std::vector<bool> gamm,
                     int depth,
                     int maxdepth,
                     std::vector<std::vector<mpf_class>>& orbit_set,
                     std::vector<std::vector<bool>>& gap,
                     std::vector<std::vector<int>>& branch_set,
                     std::vector<std::vector<bool>>& gamma_set)
{
	if (maxdepth <= depth)
	{
		orbit_set.push_back(orbit);
		gap.push_back(greedy);
		branch_set.push_back(branch_points);
		gamma_set.push_back(gamm);
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

				branch_points.push_back(i);
				std::vector<int> lobran = branch_points;

				// Get the alternate expansion.
				greedy_expand(y, beta, i+1, nbits, lorbit, gapper);

				std::vector<bool> logam = gamm;
				gamm.push_back(1);
				logam.push_back(0);

				// Recurse
				beta_expand_rec(y, beta, em, i+1, nbits, orbit, greedy,
				                branch_points, gamm,
				                depth+1, maxdepth, orbit_set,
				                gap, branch_set, gamma_set);
				beta_expand_rec(y, beta, em, i+1, nbits, lorbit, gapper,
				                lobran, logam,
				                depth+1, maxdepth, orbit_set,
				                gap, branch_set, gamma_set);
				return;
			}
		}
	}
}

void print_bits(const std::vector<bool>& bits)
{
   printf("len=%lu ", bits.size());
   for (size_t j=0; j<bits.size(); j++)
      printf("%d", (int) bits[j]);
   printf("\n");
}

void print_branches(const std::vector<int>& branches)
{
   printf("branch points len=%lu>> ", branches.size());
   for (size_t j=0; j<branches.size(); j++)
      printf("%d ", (int) branches[j]);
   printf("\n");
}

void beta_expand(mpf_class y, mpf_class beta, int em,
                 int maxdepth,
                 std::vector<std::vector<mpf_class>>& orbit_set,
                 std::vector<std::vector<bool>>& gap,
                 std::vector<std::vector<int>>& branch_set,
                 std::vector<std::vector<bool>>& gamma_set,
                 int nbits)
{
	std::vector<bool> bitseq;
	std::vector<mpf_class> orbit;
	std::vector<int> branch_points;
	std::vector<bool> gamm;
	bitseq.resize(nbits);
	orbit.resize(nbits);

	for (int i=0; i< nbits; i++)
	{
		orbit[i] = 0;
	}

	// First, get the baseline (greedy) orbit.
	greedy_expand(y, beta, 0, nbits, orbit, bitseq);

	// And now recursively get the rest.
	beta_expand_rec(y, beta, em, 0, nbits,
	                orbit, bitseq, branch_points, gamm,
	                0 /* depth*/, maxdepth,
	                orbit_set, gap, branch_set, gamma_set);

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

// Just the ordinary greedy beta expansion.
void beta_sequence(mpf_class y, mpf_class beta, int em,
                   std::vector<mpf_class>& orbit,
                   std::vector<bool>& bitseq,
                   int nbits)
{
	bitseq.resize(nbits);
	orbit.resize(nbits);

	for (int i=0; i< nbits; i++)
		orbit[i] = 0;

	// Get the baseline (greedy) orbit.
	greedy_expand(y, beta, 0, nbits, orbit, bitseq);
}

// ================================================================

// Sum the bit-sequence, returning the sum.
// This returns 0.5 * sum_i=0 b[i] (2J)^-i
//
double beta_sum(std::vector<bool> bits, double Jay)
{
	double acc = 1.0e-30;
	int bitlen = bits.size();
	if (60 < bitlen) bitlen = 60;
	for (int i=0; i<bitlen; i++)
	{
		acc *= 1.0 / (2.0*Jay);
		if (bits[bitlen-i-1])
		{
			acc += 0.5;
		}
	}
	return acc;
}

// ================================================================
