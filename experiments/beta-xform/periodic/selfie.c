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

/* Return length of bitstring. Same as ceil(log2(bitstr)). */
int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

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
	printf (" \\\\ %d", len);
	printf("%s", suf);
}

/* ================================================================= */

/* Implement the n'th golden polynomial. Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the bit-shift construction: the lowest powers of x are given by
 * right-most bits; highest powers are the left-most bits.
 */
double golden_poly(unsigned long idx, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	unsigned long bitstr = 2*idx+1;
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
static double find_zero(unsigned long idx, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = golden_poly(idx, mid);
	if (0.0 < fmid) return find_zero(idx, lo, mid);
	return find_zero(idx, mid, hi);
}

/* Return the beta value corresponding to the n'th golden polynomial.
 * Polynomial is constructed from the bit string of (2n+1).
 */
double golden_beta(unsigned long idx)
{
	return find_zero(idx, 1.0, 2.0);
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
		if (0.5 <= mid)
		{
			bitseq |= 1UL << i;

#define MIDEPSI 1.0e-15
			// Apply rounding pressure, so as to favor finite iterates
			// over periodic ones. MIDEPSI=1e-15 seems to work at low
			// orders e.g. beta from index=6.
			mid -= 0.5-MIDEPSI;
		}
		mid *= beta;
	}
	return bitseq;
}

/* Return true if `idx` is a self-describing index. This means that
 * the root of the golden_beta polynomial, when iterated by mid-point
 * iteration, reproduces the a dyadic bitstring, that, when bit-reversed,
 * gives back the index.
 *
 * The length of the orbit will be log_2(2n+1).
 *
 * This is returns exactly the same values as the theta() function,
 * except that it does not use a recursive algo to do it's work.
 */
bool valid_gold_index(unsigned long idx)
{
	double gold = golden_beta(idx);
	unsigned long dyad = beta_to_dyadic(gold);
	// printf("Index=%ld gold=%20.16g\n", idx, gold);
	// print_dyadic(dyad, 40, "Dyadic orbit=", "\n");

	unsigned long tno = 2*idx+1;
	int len = bitlen(tno);
	for (int i=0; i< len; i++)
	{
		if ((tno &1UL) != ((dyad>>(len-i-1)) &1UL)) return false;
		tno >>= 1;
	}
	return true;
}

/* ================================================================= */
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
/*
 * Given a rational p/q, return the prefix-integer, the cycle integer,
 * and the cycle length for the infinite ultimately-periodic dyadic
 * string that would result from the dyadic expansion of p/q.
 */
void get_cycle (unsigned long p, unsigned long q,
                unsigned long *pfxp, unsigned long *cycp,
                int *cyclenp)
{
	unsigned long pfx = p;
	while (pfx%2 == 0) pfx >>= 1;
	int cyclen = bitlen(q);
	unsigned long cyc = (1UL << cyclen) - q;
	*pfxp = pfx;
	*cycp = cyc;
	*cyclenp = cyclen;
}

/* --------------------------- END OF LIFE ------------------------- */
