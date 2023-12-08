/*
 * irred-fraction.c
 *
 * Find integer sequence for the golden polynomials.
 * Relate it to the continued-fraction representation.
 * See irred.c for additional utilities.
 *
 * February 2018, October 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/* Return the n'th golden polynomial. It can be constructed from
 * the bit string of (2n+1).
 */
double beta(unsigned long n, double x)
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
double find_zero(unsigned long n, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = beta(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/* Return length of bitstr, length in bits */
int len(unsigned long n)
{
	int len=0;
	unsigned long bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/** Helper array, needed for finding gold midpoints */
double* zero = NULL;
void setup_gold(int nmax)
{
	zero = (double*) malloc(nmax*sizeof(double));
	for (int i=0; i< nmax; i++) zero[i] = -1.0;
}

// Return true if this is a valid polynomial
bool zero_is_bracketed(int n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15
	// printf("---------\ngold=%g\n", gold);
	bool ork = true;
	int nhl = n;
	int nh = nhl >> 1;
	while (nh)
	{
		// printf("duuude n=%d nhl=%d nh=%d znh=%g go=%g comp=%d bad=%d\n",
		// n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}

	return ork;
}

/**
 * Find the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not really correspond to a golden
 * polynomial, return zero.
 */
double find_gold(long n)
{
	if (-1 == n) return 2.0;

	double gold = zero[n];
	if (gold < -0.5)
	{
		gold = find_zero(n, 1.0, 2.0);
		zero[n] = gold;
	}

	if (zero_is_bracketed(n, gold)) return gold;
	return 0.0;
}

// Return the real value of the corresponding continued fraction.
// This is 1/(1+m0+ 1/(1+m1+ 1/(1+m2+ ... 1/(1+mk))))
// where we have to add one because the sequence has zeros in it.
double cf_to_real(int cfrac[], int len)
{
	double ex = 0.0;
	for (int i=len-1; i >=0; i--)
	{
		ex += cfrac[i] + 1.0;
		ex = 1.0 / ex;
	}
	return ex;
}

/*
 * Given the continued-fraction representation, return the
 * corresponding sequence number.
 */
long sequence_from_cf(int cfrac[], int len)
{
	if (1 == len)
	{
		if (-1 == cfrac[0]) return -1;
		return 1UL << cfrac[0];
	}

	long leader = sequence_from_cf(cfrac, len-1);
	if (-1 == leader) return -1; // avoid overflow

	long follower = 2*leader + 1;

	// This appears to be correct .... for now ...  !?
	int shift = cfrac[len-1];

	int bump = len-3;
	if (bump < 0) bump = 0;
	for (int j=0; j<= bump; j++)
		shift += cfrac[j];

	// More careful overflow check
	int nbits = 0;
	long fo = follower;
	while (fo >>= 1) ++nbits;

	if (60 < nbits+shift) return -1;

	follower *= 1UL << shift;

	return follower;
}

void print_seq(int cfrac[], int len, char* head, char* tail)
{
	printf("%s [", head);
	for (int i=0; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

#define SZ 20
/*
 * Validate the bounds on the continued fraction.
 * ... and print it out.
 * maxn == cutoff for highest known n; this avoid overflow.
 */
void validate_cf(int cfrac[], int len, long maxn)
{
	static double prevgold = 2.0;

	long seq = sequence_from_cf(cfrac, len);
	if (seq >= maxn) return;
	if (-1 == seq) return;
	double gold = find_gold(seq);
	// printf("seq = %ld gold=%g ", seq, gold);

	double ex = cf_to_real(cfrac, len);
	printf("%ld	%g	%g	#", seq, gold, ex);
	print_seq(cfrac, len, "", "");

	if (gold >= prevgold)
		printf("FAIL total order!\n");
	else
		prevgold = gold;

	// Validate bracketing.
	if (1 < len)
	{
		print_seq(cfrac, len-1, " left=", "");

		long left = sequence_from_cf(cfrac, len-1);

		// Right bracket is tricky.
		int rfrac[SZ];
		for (int i=0; i<len; i++) rfrac[i] = cfrac[i];
		int rlen = len;
		while (1 < rlen && 0 == rfrac[rlen-1]) rlen--;

		rfrac[rlen-1]--;
		print_seq(rfrac, rlen, " right=", "");

		long right = sequence_from_cf(rfrac, rlen);

		// When bracketed by leader peers, left must be greater than right.
		if (len-1 == rlen && left <= right)
			printf("FAIL peer bracket order! left=%ld right=%ld\n", left, right);

		if (len <= rlen && 0 == rfrac[rlen-1] && right <= left && -1 != right)
			printf("FAIL hi bracket order! left=%ld right=%ld\n", left, right);

		// This check passes.... I think it's correct ...
		if (len > rlen && left <= right)
			printf("FAIL lo bracket order! left=%ld right=%ld\n", left, right);

		if (seq <= left)
			printf("FAIL left bracket! left=%ld seq=%ld\n", left, seq);

		if (seq <= right)
			printf("FAIL right bracket! right=%ld seq=%ld\n", right, seq);

		double lg = find_gold(left);
		double rg = find_gold(right);

		if (gold < lg) printf("fail left gold! left=%g gold=%g\n", lg, gold);
		if (rg < gold) printf("fail right gold! gold=%g right=%g\n", gold, rg);
	}

	printf("\n");
	fflush(stdout);
}

/*
 * Iterate on the continued fraction.
 * i.e. generate sequences
 * maxdepth == number of doubling steps
 * maxlength == max length of fraction.
 * maxn == cutoff for highest known n
 */
void iterate_cf(int cfrac[], int len, int maxdepth, int maxlength, long maxn)
{

	// Iterate to max length, first.
	if (len < maxlength)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len] = 0;
		iterate_cf(bfrac, len+1, maxdepth, maxlength, maxn);
	}

	validate_cf(cfrac, len, maxn);

	// Iterate depthwise second.
	if (cfrac[len-1] < maxdepth)
	{

		int afrac[SZ];
		for (int i=0; i<len; i++) afrac[i] = cfrac[i];
		afrac[len-1] ++;
		iterate_cf(afrac, len, maxdepth, maxlength, maxn);
	}
}

// Attempt to generate sequence expansions from integer labels.
int reverso(int cfrac[], int nseq)
{
	// int norig = nseq;

	int msum = 0;
	while (0 == nseq %2) { msum ++; nseq /=2; }

	// printf("entry norig=%d msum=%d nseq=%d\n", norig, msum, nseq);

	// Terminate recursion
	if (1 == nseq)
	{
		// printf("the end\n");
		cfrac[0] = msum;
		return 0;
	}

	// If we are here, nseq is odd; reduce it and try again.
	nseq -=1;
	nseq /=2;

	// Reject pure powers of two, when they occur after an odd number.
	if (1 < nseq && 0 == msum)
	{
		int pure = nseq;
		while (0 == pure %2) { pure /=2; }
		if (1 == pure)
		{
			// printf("ppppppppppure power!! orig=%d nseq=%d\n", norig, nseq);
			return -666;
		}
	}

	// Recurse
	int loc = reverso(cfrac, nseq);
	if (loc < 0) return loc;

	// Remove contributions from shorter sequences
	for (int j=0; j<loc; j++)
		msum -= cfrac[j];

	// Special-case the length=2 case.
	if (0 == loc) msum -= cfrac[0];

	// printf("post recur norig=%d nseq=%d loc=%d msum=%d\n", norig, nseq, loc, msum);
	cfrac[loc+1] = msum;

	return loc+1;
}

int main(int argc, char* argv[])
{
	int cfrac[SZ];

	// Attempted reverse listing.
	int nmax = 128;
	for (int n=1; n<nmax; n ++)
	{
		for (int i=0; i<SZ; i++) cfrac[i] = -666;
		int len = 1 + reverso(cfrac, n);
		if (len < 0)
		{
			// printf(">>>>> %d rejected\n", n);
			continue;
		}
		int seqno = sequence_from_cf(cfrac, len);
		if (n != seqno)
		{
			printf("Sequence numbering fail!! in=%d out=%d", n, seqno);
		}
		printf(">>>>> %d len=%d ", n, len);
		print_seq(cfrac, len, "sequence ", "\n");
	}

// #define SANITY_CHECK
#ifdef SANITY_CHECK
	// Sanity check, only; short runtime.
	int norder = 18;
	int nmax = (1<<norder) + 1;

	setup_gold(nmax);

	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	iterate_cf(cfrac, 1, 3, 8, nmax);
#endif

// #define BIG_GRAPH
#ifdef BIG_GRAPH
	// Using 1<<24 takes about 50 seconds to setup gold.
	int norder = 26;
	int nmax = (1<<norder) + 1;

	setup_gold(nmax);

	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	// iterate_cf(cfrac, 1, 10, 10, nmax);

	// Need to go to high depth to avoid big gap at golden mean.
	int maxdepth = 16;
	int maxlen = 9;
	printf("#\n# Max order of polynomials = %d num=2^order = %d\n", norder, nmax);
	printf("#\n# Iterate to maxdepth=%d maxlen=%d\n#\n", maxdepth, maxlen);
	fflush (stdout);
	iterate_cf(cfrac, 1, maxdepth, maxlen, nmax);
#endif
}
