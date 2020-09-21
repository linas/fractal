/*
 * sidorov.C
 *
 * Understand bit-sequence mappings.  The expander and compresor
 * functions. This time with alternative expansions, per sidorov.
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#include <vector>

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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

#define NBITS 50

// Compute two alternate beta expansions, and compare them
// side by side. This is currently for debugging only; it
// spews too many prints. Use beta_expand below for the real
// thing.
//
// Note that beta = 2*K
double sdr(double y, double K, int em)
{
	// Generate the beta expansion bits in a greedy fashion.
	char grebits[NBITS];
	double greedy[NBITS];
	for (int i=0; i<NBITS; i++)
	{
		greedy[i] = y;
		if (0.5 <= y)
		{
			y -= 0.5;
			grebits[i] = 1;
		}
		else grebits[i] = 0;
		y *= 2.0*K;
	}

	// Search for em runs.
	// The Sidorov paper has an error, the m is off by one.
	char lobits[NBITS];
	for (int i=0; i<NBITS-em; i++)
	{
		lobits[i] = grebits[i];
		if (1 == grebits[i])
		{
			bool found = true;
			for (int j=1; j<=em; j++)
			{
				if (1 == grebits[i+j])
				{
					found = false;
					break;
				}
			}
			if (found)
			{
				lobits[i] = 0;
				y = greedy[i];
				printf("# got one at i=%d y=%g next=%g\n", i, y, y*2*K);
				y *= 2.0*K;

				for (int j=i+1; j<NBITS; j++)
				{
					if (0.5 <= y)
					{
						y -= 0.5;
						lobits[j] = 1;
					}
					else lobits[j] = 0;
					y *= 2.0*K;
				}
				break;
			}
		}
	}

	double Jay = K;

	// Reconstruct both sequences;
	double hiacc = 1.0e-30;
	double loacc = 1.0e-30;
	for (int i=0; i<NBITS; i++)
	{
		hiacc *= 1.0 / (2.0*Jay);
		loacc *= 1.0 / (2.0*Jay);
		if (grebits[NBITS-i-1])
		{
			hiacc += 0.5;
		}
		if (lobits[NBITS-i-1])
		{
			loacc += 0.5;
		}
	}

	printf("# hi=");
	for (int i=0; i<NBITS; i++) printf("%d", grebits[i]);
	printf("\n");
	printf("# lo=");
	for (int i=0; i<NBITS; i++) printf("%d", lobits[i]);
	printf("\n");

	return hiacc-loacc;
}

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

// Generate a collection of equivalent beta expansions at K.
// `y` is the number to expand.
// `K` is beta/2
// `em` is the number of zero bits to look for.
std::vector<std::vector<bool>> beta_expand(double y, double K, int em)
{
	std::vector<std::vector<bool>> bitset;

	// First, get the baseline orbit.
	std::vector<double> orbit;
	std::vector<bool> greedy;
	orbit.resize(NBITS);
	greedy.resize(NBITS);

	greedy_expand(y, K, 0, orbit, greedy);
	bitset.push_back(greedy);

	// Search for em runs.
	// The Sidorov paper has an error, the m is off by one.
	for (int i=0; i<NBITS-em; i++)
	{
		std::vector<bool> gapper = greedy;
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
				// Set to zero, and resume expansion.
				gapper[i] = 0;
				y = orbit[i];
				y *= 2.0*K;
				orbit[i] = y;

				// Get the alternate expansion.
				greedy_expand(y, K, i+1, orbit, gapper);
				bitset.push_back(gapper);
				greedy = gapper;
			}
		}
	}

	// Debug print
	for (size_t n=0; n<bitset.size(); n++)
	{
		std::vector<bool> bits = bitset[n];
		printf("# n=%lu ", n);
		for (int i=0; i<NBITS; i++)
			printf("%d", (int) bits[i]);
		printf("\n");
	}
	return bitset;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d\n#\n", Kay, em);

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

	int npts = 13;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		beta_expand(x, Kay, em);
		printf("%d	%g\n", i, x);
	}
}
