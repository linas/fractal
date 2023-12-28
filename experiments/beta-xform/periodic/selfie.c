/*
 * selfie.c
 *
 * Redesign of polynomial encoding and root finding. The primary benefit
 * of this redesign is that intermediate beta values are not cached, and
 * so this does not consume RAM. Nor is it any slower for most tasks, so
 * it's a win over-all.
 *
 * Includes code for the finite orbits, only. See selfie-rational for the
 * ultimately-periodic orbits.
 *
 * This is a redesign and sometimes copy of code in irred-gold.c
 *
 * Notes about bit-encodings.
 * --------------------------
 * Four different bit-encodings are in use. They are:
 *
 * -- "Index encoding". A long unsigned int.  Used for finite orbits.
 *    Not all indexes are valid. Use valid_gold_index() to determine
 *    if it is OK.
 *
 * -- "Leading-one scheme" aka "big-endian encoding", aka "finite orbit
 *    encoding". Used for finite orbits. MSB is to far left, LSB is far
 *    right, both MSB and LSB are one.  No length specifier needed. MSB
 *    is always one, because first move of orbit is always one. LSB is
 *    always one, because orbit always terminates. Equal to (2*idx+1).
 *
 * -- "Dyadic encoding" aka "2-adic encoding" aka "little-endian"
 *    encoding. Used for unbounded orbits. MSB is far right, later
 *    bits in orbit are to the left. Explicit length is needed.
 *    Midpoint orbits always begin with a 1, so the right-most bit
 *    is always a one. Bit-reversed version of big-endian encoding,
 *    with length suaitably indicated.
 *
 * -- "Tree encoding", used to indicate location in a full binary tree.
 *    First move is the right-most bit, and it must always be a one.
 *    Later moves down the tree are bits to the left of that. Explicit
 *    length is needed. This is more-or-less the same thing as the
 *    2-adic encoding. Just a different interpretation.
 *
 *    This encoding provides the canonical index numbering for a tree.
 *    For every non-negative integer, there is a unique tree location.
 *    0 == before the root, to the left
 *    1 == root node
 *    2,3 == first row below the root, L and R
 *    4,5,6,7 == next row down.
 *
 * The big-endian encoding is used for finite-length orbits. It places
 * the MSB to far-left. The left-most bit is always one. Usually, so is the
 * right-most bit. Some comments call this the leading-one scheme.
 *
 * Advantages:
 * -- No length specification is needed. The leading-one indicates the
 *    start of the string.
 * -- Can be interpreted as an integer: thus, finite orbits can be
 *    associated with an explicit integer.
 *
 * Disadvantages:
 * -- Long orbits have extremely large integers. Periodic orbits cannot
 *    be represented at all.
 * -- "nearby" orbits, sharing common leading bitstrings, have completely
 *    different integer values.
 *
 * The little-endian encoding is used for dyadic fractions and for
 * orbits of unbounded length. This is the 2-adic encoding; infinite
 * strings extend to the left.  The MSB is to the far-right, so the
 * MSB of the dyadic is the LSB of C/C++. (This is normal for p-adics.)
 *
 * The 2-adic encoding means:
 *   -- bit-zero is the right-most bit.
 *   -- bit-one is <<1 to the left of bit zero.
 *   -- finite-length strings need an explicit length.
 *   -- printing will be left-to-right, so only the storage of
 *      the string is 2-adic, but not the printouts.
 *
 * The two different encodings can cause some confusion. As a general
 * rule: if the orbit is finite, use the leading-one scheme, and get an
 * integer index for that orbit. If the orbit is periodic or infinite,
 * use the 2-adic encoding, as it seems easier to deal with, although
 * it does require passing around an explicit length, when a length is
 * needed. Phew. Glad we got that straight. I was getting worried there,
 * for a while.
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define WORDLEN 64
#define MAXIDX (1UL<<(WORDLEN-1))
#define NEG_ONE ((unsigned long)(-1L))

/* Return length of bitstring. Same as ceil(log2(bitstr)). */
int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

// Print bitstring, so that bit zero appears left-most.
// Roughly speaking, this reverses the bit-order.
void print_dyadic(unsigned long bitseq, int len,
                  const char* pre, const char* suf)
{
	printf("%s", pre);
	int i;
	for (i=0; i<len; i++)
	{
		printf("%ld", bitseq & 1UL);
		bitseq >>= 1;
		if (0 == bitseq) break;
	}
	printf (" \\\\ %d", i+1);
	printf("%s", suf);
}

/*
 * Return dyadic string corresponding to beta. This is obtained by
 * midpoint iteration, storing the bits in dyadic order, so that
 * first iteration is right-most bit. Note that right-most bit will
 * always be one, since the mid-point is always greater than 1/2.
 *
 * Finite orbits will return a finite string. Chaotic orbits might
 * use all 64 bits in the unsigned long.
 *
 * The bitstring that is returned here is an the valid integer beta-
 * index string for (2n+1) but in reversed order.
 */
unsigned long beta_to_dyadic(double beta)
{
	unsigned long bitseq = 0;
	double mid = 0.5*beta;
	int numzero=0;
	for (int i=0; i < WORDLEN; i++)
	{
		if (0.5 <= mid)
		{
			numzero = 0;
			bitseq |= 1UL << i;

#define MIDEPSI 1.0e-15
			// Apply rounding pressure, so as to favor finite iterates
			// over periodic ones. Except this doesn't work as expected.
			// MIDEPSI=1e-15 gets 3 to work, but then 6 doesn't. And other
			// very mixed, unpredictable results. Bummer. However, this is
			// delicate, and if we don't use MIDEPSI, then the self-description
			// mechanism fails. So don't just tinker with this without running
			// unit tests.
			mid -= 0.5-MIDEPSI;
			// mid -= 0.5;
		}
		else
		{
			numzero++;
			// If we have a string of 50 zeros, assume a rounding error,
			// and bail out. This avoids some of the reconstruction failures.
			// For example, index 3 fails at bit 54
			// index 13 fails at bit 61
			if (WORDLEN-15 < numzero) break;
		}
		mid *= beta;
	}
	return bitseq;
}

// Just as above, but reversing the bit-sequence, and then
// dividing by two to get the index.
unsigned long beta_to_index(double beta)
{
	unsigned long orbit = beta_to_dyadic(beta);
	unsigned long idx = 0UL;
	for (int i=0; i<WORDLEN; i++)
	{
		idx <<= 1;
		idx |= orbit & 1UL;
		orbit >>= 1;
		if (0 == orbit) break;
	}
	idx >>= 1;
	return idx;
}

/* ================================================================= */

/* Implement the n'th golden polynomial. Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the bit-shift construction: the lowest powers of x are given by
 * right-most bits; highest powers are the left-most bits.
 *
 * Should return exactly the same values as the recursive version below.
 */
double xgolden_poly(unsigned long idx, double x)
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
	// printf("bisect n=%ld x=%18.16g beta=%g\n", idx, x, xn-acc);
	return xn - acc;
}

/* Implement the n'th golden polynomial. Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the recursive construction: recurse according to the 2-adic bit position.
 * Should provide exactly the same results as above.
 */
double golden_poly(unsigned long idx, double x)
{
	double p_n = x-1.0;   // start with p_0(x)
	unsigned long bitstr = idx;
	int len = bitlen(bitstr);
	for (int i=0; i< len; i++)
	{
		if (bitstr & (1UL<<(len-1-i)))
			p_n = x * p_n - 1.0;
		else
			p_n = x * (p_n + 1.0) - 1.0;
	}
	// printf("bisect n=%ld x=%18.16g beta=%g\n", idx, x, p_n);
	return p_n;
}

/* Use midpoint bisection to find the single, unique
 * positive real zero of the n'th golden polynomial.
 */
static double find_zero(unsigned long idx, double lo, double hi)
{
	int nextra = 0;
	while (1)
	{
		double mid = 0.5 * (lo+hi);
		if (1.0e-15 > hi-lo) nextra++;
		if (3 < nextra) return mid; // run the loop three more times.
		double fmid = golden_poly(idx, mid);
		if (0.0 == fmid) return mid;
		if (0.0 < fmid)
			hi = mid;
		else
			lo = mid;
	}
}

/* Return the beta value corresponding to the n'th golden polynomial.
 * Polynomial is constructed from the bit string of (2n+1).
 */
double golden_beta(unsigned long idx)
{
	if (NEG_ONE == idx) return 1.0;
	if (0 == idx) return 2.0;
	return find_zero(idx, 1.0, 2.0);
}

/* ================================================================= */
/* Return true if `idx` is a self-describing index. This means that
 * the root of the golden_beta polynomial, when iterated by mid-point
 * iteration, reproduces the a dyadic bitstring, that, when bit-reversed,
 * gives back the index.
 *
 * The length of the orbit will be log_2(2n+1).
 *
 * This is returns exactly the same values as the theta() function,
 * except that it does not use a recursive algo to do it's work.
 * Modernized version of is_valid_index().
 *
 * Caution: This is somewhat deleicate w.r.t. the midpoint iteration.
 * rounding errors can result in a midpoint that fails to terminate
 * within rounding errors, which results in nearby orbits that give the
 * correct beta, but fail the exact bit-pattern match, below. Hmmm.
 */
bool valid_gold_index(unsigned long idx)
{
	// Far left limit is beta=1 encoded as -1 and it is valid.
	if (NEG_ONE == idx) return true;

	double gold = golden_beta(idx);
	unsigned long dyad = beta_to_dyadic(gold);
	// printf("Index=%ld gold=%20.16g\n", idx, gold);
	// print_dyadic(dyad, 40, "Dyadic orbit=", "\n");

	// Verify that the reconstructed bit-sequence reproduces the index
	// we are given.
	unsigned long tno = 2UL * idx + 1UL;
	int len = bitlen(tno);
	for (int i=0; i< len; i++)
	{
		if ((tno & 1UL) != ((dyad>>(len-i-1)) & 1UL)) return false;
		tno >>= 1;
	}
	return true;
}

// Find the leader of a valid index.  That is, find the index of the form
// 2^h(2k+1) with the smallest height h that gives a valid index.
// This move can be undone with bracket_gold_left()
// This is the modernized version of find_leader()
unsigned long gold_leader(unsigned long idx)
{
	// Extreme left side is marked with -1. So jump to center, which is 1.
	if (NEG_ONE == idx) return 1;

	// Indicate overflow errors
	if (MAXIDX <= idx) return -2L;

	idx = 2UL * idx + 1UL;
	while (false == valid_gold_index(idx))
	{
		idx <<= 1;
		if (MAXIDX <= idx) return -2L;  // overflow error
	}
	return idx;
}

// Same as above; make it clear it is a right-move.
// This move can be undone with bracket_gold_left()
unsigned long move_gold_right(unsigned long idx)
{
	return gold_leader(idx);
}

// Move left on the binary tree. This is easy!
// This move can be undone with bracket_gold_right()
unsigned long move_gold_left(unsigned long idx)
{
	// Extreme right side is zero. Move to 1 in the center.
	if (0 == idx) return 1UL;

	if (MAXIDX <= idx) return -2L;
	return idx << 1;
}

// Given an valid index, return the right side of the bracket that
// contains that index. This undoes what move_gold_right() did.
unsigned long bracket_gold_left(unsigned long idx)
{
	unsigned long clef = idx;
	while (0 == clef%2 && 0 != clef) clef >>= 1;
	clef = (clef - 1UL) / 2UL;
	if (0 == clef) clef = (unsigned long) -1L; // Yes really.
	return clef;
}

// Given an valid index, return the right side of the bracket that
// contains that index. This undoes what move_gold_left() did.
unsigned long bracket_gold_right(unsigned long idx)
{
	unsigned long brig = idx >> 1;

	// If it's odd, recurse.
	if (1 == idx%2)
		return bracket_gold_right(brig);

	// If its good, we are done.
	if (valid_gold_index(brig)) return brig;

	// Ooops. Must have been a leader. Recurse.
	while (brig%2 == 0) brig >>= 1;
	return bracket_gold_right(brig);
}

/* ================================================================= */
/*
 * Simple, basic unit-test for indexes.
 */
bool test_gold_index(unsigned long idx)
{
	if (MAXIDX <= idx)
		{ printf("Error: overflow index %ld\n", idx); return false; }

	bool ok = valid_gold_index(idx);
	if (!ok)
		printf("Error: Not a valid index: %ld\n", idx);

	double beta = golden_beta(idx);

	// ----------------------
	long cleft = bracket_gold_left(idx);
	bool leftok = valid_gold_index(cleft);
	if (!leftok)
		{ printf("Error: Left bracket %ld of %ld not valid\n", cleft, idx); ok=false; }

	double gleft = golden_beta(cleft);
	if (gleft >= beta)
		{ printf("Error: bad left bracket at %ld: nleft=%ld gold=%g gleft=%g\n",
			idx, cleft, beta, gleft); ok = false; }

	// ----------------------
	long cright = bracket_gold_right(idx);
	bool rightok = valid_gold_index(cright);
	if (!rightok)
		{ printf("Error: Left bracket not valid: %ld\n", idx); ok=false; }

	double gright = golden_beta(cright);
	if (gright < beta)
		{ printf("Error: bad right bracket at %ld: nright=%ld gold=%g gright=%g\n",
			idx, cright, beta, gright); ok = false; }

	// ----------------------
	// Validate dyadic moves.
	unsigned long dleft = move_gold_left(idx);
	unsigned long rup = bracket_gold_right(dleft);
	if (rup != idx)
		{ printf("Error: Failed return from left move for %ld went to %ld came back to %ld\n",
			idx, dleft, rup); ok = false; }

	unsigned long dright = move_gold_right(idx);
	unsigned long lup = bracket_gold_left(dright);
	if (lup != idx)
		{ printf("Error: Failed return from right move for %ld went to %ld came back to %ld\n",
			idx, dright, lup); ok = false; }

	if (!ok) printf("-------\n");
	return ok;
}

/* Print some information about this index */
void print_gold_info(long seqno)
{
	double beta = golden_beta(seqno);
	printf("Index %ld beta = %18.16g\n", seqno, beta);
	unsigned long dyad = beta_to_dyadic(beta);
	print_dyadic(dyad, 40, "midpo orbi= ", "\n");

	test_gold_index(seqno);

	long nleft = bracket_gold_left(seqno);
	double gleft = golden_beta(nleft);
	printf("Left limit: %ld = %g\n", nleft, gleft);
	unsigned long ldyad = beta_to_dyadic(gleft);
	print_dyadic(ldyad, 40, "left  orbi= ", "\n");

	long nright = bracket_gold_right(seqno);
	double gright = golden_beta(nright);
	printf("Right limit: %ld = %g\n", nright, gright);
	unsigned long rdyad = beta_to_dyadic(gright);
	print_dyadic(rdyad, 40, "right orbi= ", "\n");
	printf("-------\n");
}

/* --------------------------- END OF LIFE ------------------------- */
