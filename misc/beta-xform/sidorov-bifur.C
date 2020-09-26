/*
 * sidorov-density.C
 * Graph fake bifurcation diagrams.
 *
 * Linas Vepstas Sept 2020
 */

#include "brat.h"

#include "sidorov-big.C"


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
#define MAXDEPTH 12

	static bool init=false;
	if (not init)
	{
		do_init(nbits);
		init = true;

double Kay=0.83;
double x=param;
		make_random_bitsequence(beta, 2.0*Kay, nbits, 1.0e9);
		int em = emrun(Kay);

		mpf_class ex;
		make_random_bitsequence(ex, x, nbits, 1e9);

		beta_expand(ex, beta, em, MAXDEPTH, orbit_set, bitset, branch_set, nbits);
	}

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double Jay = 1.0;
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
			if (bitseq[bitlen-i-1])
			{
				acc += 0.5;
			}
		}

		acc *= array_size;
		int pix = acc;
		if (array_size <= pix) pix = array_size-1;

		array[pix] = 1.0;
	}
}

DECL_MAKE_BIFUR(fake_bifur)

// ================================================================
