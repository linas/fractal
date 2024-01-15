/*
 * deriv.c
 *
 * Explorederivative of polynomials.
 *
 * January 2024
 */

#include "selfie.c"

/* ================================================================= */

/* Implement derivative of the n'th golden polynomial.
 * Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the bit-shift construction: the lowest powers of x are given by
 * right-most bits; highest powers are the left-most bits.
 */
double golden_poly_deriv(unsigned long idx, double x)
{
	double acc = 0.0;
	double xn = 1.0 / x;
	double n = 0.0;
	unsigned long bitstr = 2*idx+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += n * xn;
		xn *= x;
		n += 1.0;
		bitstr >>= 1;
	}

	return n*xn - acc;
}

/* ================================================================= */

int main(int argc, char* argv[])
{
#define MANUAL_EXPLORER
#ifdef MANUAL_EXPLORER
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxorder>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	long maxidx = 1UL << maxord;

	for (long idx = 1; idx < maxidx; idx++)
	{
		bool ok = valid_gold_index(idx);
		if (!ok) continue;

		int ord = bitlen(idx);
		long twoo = 1UL << ord;
		double gold = golden_beta(idx);
		double deriv = golden_poly_deriv(idx, gold);
		deriv /= twoo;
		// deriv = 1.0 / deriv;
		printf("%ld	%g	%g\n", idx, gold, deriv);
	}
#endif
}

/* --------------------------- END OF LIFE ------------------------- */
