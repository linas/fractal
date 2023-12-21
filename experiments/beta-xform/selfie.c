/*
 * selfie.c
 *
 * Redesign of polynomial encoding and root finding.
 *
 * Notes about bit-encodings.
 * Orbits and dyadic fractions will use a little-endian encoding.
 * This is the same as 2-adic encoding, so that infinite strings
 * extend to the left.
 *
 * This means:
 *   -- bit-zero is the right-most bit.
 *   -- bit-one is <<1 to the left of bit zero.
 *   -- finite-length strings need an explicit length.
 *   -- printing will be left-to-right, so only the storage of
 *      the string is 2-adic, but not the printouts.
 *
 * This is a redesign and sometimes copy of code in irred-gold.c
 *
 * December 2023
 */

// Print bitstring, so that bit zero appears left-most.
// Roughly speaking, this reverses the bit-order.
void print_dyadic(unsigned long bitseq, int len, char* pre, char* suf)
{
	printf("%s", pre);
	for (int i=0; i<len; i++)
	{
		printf("%ld", bitseq & 1UL);
		bitseq >>= 1;
	}
	printf (" \\%d", len);
	printf("%s", suf);
}

// Return the dyadic that approximates the rational p/q,
// truncated to len bits. Rationals have infinite periodic bitseqs,
// so this just expands that bitseq, and then truncates it.
//
// Bit strings are stored right-to-left, so that the decimal point at
// far right. This is a 2-adic storage convention. For example, 1/2 and
// 1/4 and 1/8 are all stored as 0x1 but have length 1,2,3.
//
unsigned long rational_to_dyadic(unsigned long p, unsigned long q, int len)
{
	unsigned long bitseq = 0;
	for (int i=0; i<len; i++)
	{
		bitseq <<= 1;
		p <<= 1;
		if (q <= p)
		{
			bitseq |= 1UL;
			p -= q;
		}
	}
	return bitseq;
}

/* ================================================================= */

/* Implement the n'th golden polynomial. Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the bit-shift construction: the lowest powers of x are given by
 * right-most bits; highest powers are the left-most bits.
 */
double golden_poly(unsigned long n, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	unsigned long bitstr = 2*n+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += xn;
		xn *= x;
		bitstr >>= 1;
	}
// printf("duuude n=%d x=%20.16g beta=\n", n, x, xn-acc);
	return xn - acc;
}

/* Use midpoint bisection to find the single, unique
 * positive real zero of the n'th golden polynomial.
 */
static double find_zero(unsigned long n, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = golden_poly(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/* Return the beta value corresponding to the n'th golden polynomial.
 * Polynomial is constructed from the bit string of (2n+1).
 */
double golden_beta(unsigned long n)
{
	return find_zero(n, 1.0, 2.0);
}

/* ================================================================= */

/*
 * Return dyadic string corresponding to beta. This is obtained by
 * midpoint iteration, storing the bits in dyadic order, so that
 * first iteration is right-most bit. Note that right-most bit will
 * always be one, since the mid-point is always greater than 1/2.
 *
 * Finite orbits will return a finite string. Chaotic orbits might
 * use all 64 bits in the unsigned long.
 *
 * "Valid" integer betas will return the bitstring for (2n+1) but
 * in reversed order.
 */
unsigned long beta_to_dyadic(double beta)
{
	unsigned long bitseq = 0;
	double mid = 0.5*beta;
	for (int i=0; i < 8*sizeof(unsigned long); i++)
	{
		if (0.5 < mid)
		{
			bitseq |= 1UL << i;
			mid -= 0.5;
		}
		mid *= beta;
	}
	return bitseq;
}

/* Construction
 * is the mid-point construction: repeated iteration of the midpoint
 * 1/2 with this beta will (re-)generate the same bitstring, until
 * returning to the midpoint. The bits are just whether the orbit went
 * left or right of midpoint. The length of the orbit will be log_2(2n+1).
 */
/* --------------------------- END OF LIFE ------------------------- */
