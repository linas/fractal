/*
 * sidorov-density.C
 * Graph fake bifurcation diagrams.
 * Great hopes for this fizzled. There are trees, but not visually dramatic.
 *
 * Linas Vepstas Sept 2020
 */

#include "brat.h"

#include "sidorov-big.C"

double wave(int k, double front)
{
	if (front < k) return 0.0;

	double hard=3.0;
	// double hard=front;
	// if (hard < front-k) return 1.0;

	// double back = (front - k) / front;  // zero to one.
	double back = (front - k) / hard;  // zero to one.

	// return sqrt(back);
	// return back;
	// return back*back;

	if (1.0< back) return 1.0;
	return back;

	// Smooth transition
	double bump = exp(-0.25/(back*back));
	back = 1.0-back;
	double flip = exp(-0.25/(back*back));

	double trans = bump / (flip+bump);
	return trans;
}

static void fake_bifur (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double param)
{
	static mpf_class beta;
	int nbits = itermax;

	static std::vector<std::vector<mpf_class>> orbit_set;
	static std::vector<std::vector<bool>> bitset;
	static std::vector<std::vector<int>> branch_set;

// #define MAXDEPTH 16
#define MAXDEPTH 7

double Kay=0.83;
// Kay = 0.55;
Kay = 0.75;
double x=param;
	static bool init=false;
	if (not init)
	{
		do_init(nbits);
		init = true;

		make_random_bitsequence(beta, 2.0*Kay, nbits, 1.0e9);
		int em = emrun(Kay);

		mpf_class ex;
		make_random_bitsequence(ex, x, nbits, 1e9);

		beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);

#ifdef PRINT_BITS
		int npaths = bitset.size();
		for (int ipath = 0; ipath<npaths; ipath++)
		{
			std::vector<bool> bitseq = bitset[ipath];
			printf("%d	", ipath);
			int bitlen = bitseq.size();
			if (60 < bitlen) bitlen = 60;
			for (int i=0; i<bitlen; i++)
				printf("%d", (int)bitseq[i]);
			printf("\n");
		}
#endif
	}

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

// printf("duude row=%g\n", row);
	double Jay = 1.0;
// Jay = 0.7;
	double front = pow(2, row);
	int npaths = bitset.size();
	for (int ipath = 0; ipath<npaths; ipath++)
	{
		std::vector<bool> bitseq = bitset[ipath];

		double acc = 1.0e-30;
		int bitlen = bitseq.size();
		if (60 < bitlen) bitlen = 60;
		for (int i=0; i<bitlen; i++)
		{
			acc *= 1.0 / (2.0*Jay);
			int kbit = bitlen-i-1;
			if (bitseq[kbit])
			{
				acc += 0.5 * wave(kbit, front);
			}
			else
			{
				// acc -= 0.1 * wave(kbit, front);
			}
		}

		acc *= array_size;
		int pix = acc;
		if (array_size <= pix) pix = array_size-1;
		if (pix < 0) pix = 0;

		array[pix] = 1.0;
	}
}

DECL_MAKE_BIFUR(fake_bifur)

// ================================================================
