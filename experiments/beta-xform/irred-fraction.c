/*
 * irred-fraction.c
 *
 * Find integer sequences for the golden polynomials.
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

// Reverse order of digits in cfrac, left to right.
void reverse_cf(int rfrac[], int cfrac[], int len)
{
	if (1 >= len) return;
	for (int i=0; i<len; i++) rfrac[i] = cfrac[len-1-i];
}

// Reverse the digits first, then get the cf value
double reverse_cf_to_real(int cfrac[], int len)
{
	int rfrac[len];
	reverse_cf(rfrac, cfrac, len);
	return cf_to_real(rfrac, len);
}

/*
 * Given a finite-lengtth Baire representation, return the
 * corresponding index number.
 */
long index_from_fbaire(int cfrac[], int len)
{
	if (1 == len)
	{
		if (-1 == cfrac[0]) return -1;
		return 1UL << cfrac[0];
	}

	long leader = index_from_fbaire(cfrac, len-1);
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

// Generate finite-Baire sequence expansions from integer index.
// Given an index, set "cfrac" to the matching sequence.
// Return the length of the cfrac sequence.
// This is the inverse of what index_from_fbaire
int index_to_fbaire(int cfrac[], unsigned long nseq)
{
	#define DBG(X) printf X
	DBG(("enter nseq=%ld\n", nseq));
	int msum = 0;
	while (0 == nseq %2) { msum ++; nseq /=2; }

	// Terminate recursion
	if (1 == nseq)
	{
		DBG(("the end msum=%d\n", msum));
		cfrac[0] = msum;
		return 1;
	}

	// If we are here, nseq is odd; reduce it and try again.
	nseq = (nseq-1)/2;

	DBG(("red after odd to nseq=%ld msum=%d\n", nseq, msum));
	// Reject pure powers of two, when they occur after an odd number.
	// Long series must terminate with an odd number at the bottom.
	if (1 < nseq)
	{
		int pure = nseq;
		while (0 == pure %2) { pure /=2; }
		if (1 == pure)
			return -666;
	}

	// Recurse.
	int len = index_to_fbaire(cfrac, nseq);

	DBG(("post recursion for nseq=%ld msum=%d new len=%d\n", nseq, msum, len));
	// Pass rejection slips back up the chain.
	if (len < 0) return len;

	// Remove contributions from shorter sequences
	for (int j=0; j<len-2; j++)
		msum -= cfrac[j];

	// Special-case the length=1 case.
	if (1 == len) msum -= cfrac[0];

	DBG(("post subtract for nseq=%ld msum=%d new len=%d\n", nseq, msum, len));
	// If more powers of two removed than can exist, then reject.
	if (msum < 0) return -555;

	// What's left over is m_k
	cfrac[len] = msum;

	// Return the length of the beast.
	return len+1;
}

void print_seq(int cfrac[], int len, char* head, char* tail)
{
	printf("%s [", head);
	for (int i=0; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

#define SZ 60
/*
 * Validate the bounds on the finite-Baire sequence representation.
 * ... and print it out.
 * maxn == cutoff for highest known n; this avoids overflow.
 */
void validate_fbaire(int cfrac[], int len, long maxn)
{
	static double prevgold = 2.0;

	long seq = index_from_fbaire(cfrac, len);
	if (seq >= maxn) return;
	if (-1 == seq) return;
	double gold = find_gold(seq);
	// printf("seq = %ld gold=%g ", seq, gold);

#define PRINT_DATA_FOR_GRAPH
#ifdef PRINT_DATA_FOR_GRAPH
	double ex = cf_to_real(cfrac, len);
	double rex = reverse_cf_to_real(cfrac, len);
	printf("%ld	%g	%g	%g #", seq, gold, ex, rex);
	print_seq(cfrac, len, "", "");
#endif

	if (gold >= prevgold)
		printf("FAIL total order!\n");
	else
		prevgold = gold;

	// Validate bracketing.
	if (1 < len)
	{
#ifdef PRINT_DATA_FOR_GRAPH
		print_seq(cfrac, len-1, " left=", "");
#endif

		long left = index_from_fbaire(cfrac, len-1);

		// Right bracket is tricky.
		int rfrac[SZ];
		for (int i=0; i<len; i++) rfrac[i] = cfrac[i];
		int rlen = len;
		while (1 < rlen && 0 == rfrac[rlen-1]) rlen--;

		rfrac[rlen-1]--;
#ifdef PRINT_DATA_FOR_GRAPH
		print_seq(rfrac, rlen, " right=", "");
#endif

		long right = index_from_fbaire(rfrac, rlen);

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

	int seqno = index_from_fbaire(cfrac, len);
	int dfrac[SZ];
	for (int i=0; i<SZ; i++) dfrac[i] = -666;
	int dlen = index_to_fbaire(dfrac, seqno);
	if (len != dlen)
		printf("Error: length violation!\n");
	for (int i=0; i<len; i++)
	{
		if (cfrac[i] != dfrac[i])
			printf("Error: sequence violation!\n");
	}

#ifdef PRINT_DATA_FOR_GRAPH
	printf("\n");
	fflush(stdout);
#endif
}

/*
 * Generate correctly-ordered finite-baire sequences. The ordering
 * is such that the corresponding golden number is strictly decreasing
 * for the generated cfrac sequences.
 *
 * maxdepth == number of doubling steps
 * maxlength == max length of fraction.
 * maxn == cutoff for highest known n
 */
void iterate_fbaire(int cfrac[], int len, int maxdepth, int maxlength, long maxn)
{
	// Iterate to max length, first.
	if (len < maxlength)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len] = 0;
		iterate_fbaire(bfrac, len+1, maxdepth, maxlength, maxn);
	}

	validate_fbaire(cfrac, len, maxn);

	// Iterate depthwise second.
	if (cfrac[len-1] < maxdepth)
	{
		int afrac[SZ];
		for (int i=0; i<len; i++) afrac[i] = cfrac[i];
		afrac[len-1] ++;
		iterate_fbaire(afrac, len, maxdepth, maxlength, maxn);
	}
}

int main(int argc, char* argv[])
{
#if 1
	int cfrac[SZ];
	for (int i=0; i<SZ; i++) cfrac[i] = -666;
	int nind = 27;
	int len = index_to_fbaire(cfrac, nind);
	printf("Index: %d len=%d", nind, len);
	print_seq(cfrac, len, " ", "\n");
#endif

// #define MANUAL_EXPLORER
#ifdef MANUAL_EXPLORER
	// Obtain one sequence from command line. Print it's index.
	if (1 == argc) {
		fprintf(stderr, "Usage: %s <sequence>\n", argv[0]);
		exit(1);
	}
	int len = argc-1;
	int cfrac[SZ];
	for (int i=0; i<len; i++) cfrac[i] = atoi(argv[i+1]);
	int seqno = index_from_fbaire(cfrac, len);
	printf("Index: %d len=%d", seqno, len);
	print_seq(cfrac, len, " ", "\n");
#endif

// #define VERFIY_FBAIRE
#ifdef VERFIY_FBAIRE

	// Verify reverse listing.
	// Takes 30 cpu-seconds to get to 1<<26
	int nmax = 1<<26;
nmax=64;
	for (int n=1; n<nmax; n ++)
	{
		int cfrac[SZ];
		for (int i=0; i<SZ; i++) cfrac[i] = -666;
		int len = index_to_fbaire(cfrac, n);
		if (len < 0)
		{
#define PRINT_SEQS
#ifdef PRINT_SEQS
			printf(">>>>> %d rejected\n", n);
#endif
			continue;
		}
		int seqno = index_from_fbaire(cfrac, len);
		if (n != seqno)
		{
			printf("Sequence numbering fail!! in=%d out=%d", n, seqno);
		}
#ifdef PRINT_SEQS
		printf(">>>>> %d len=%d ", n, len);
		print_seq(cfrac, len, "sequence ", "\n");
#endif
	}
	printf("Verified everything up to nmax=%d\n", nmax);
#endif

// #define REVERSE_VERIFY
#ifdef REVERSE_VERIFY
	// Verify that the reverse-indexing works.
	int norder = 20;
	int nmax = (1<<norder) + 1;

	setup_gold(nmax);

	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	iterate_fbaire(cfrac, 1, 5, 8, nmax);
#endif

// #define SANITY_CHECK
#ifdef SANITY_CHECK
	// Print the finite-Baire sequences. Sanity check, only; short runtime.
	int norder = 18;
	int nmax = (1<<norder) + 1;

	setup_gold(nmax);

	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	iterate_fbaire(cfrac, 1, 3, 8, nmax);
#endif

// #define BIG_GRAPH
#ifdef BIG_GRAPH
	// Obtain order from command line.
	if (4 != argc) {
		fprintf(stderr, "Usage: %s <order> <maxdepth> <maxlen>\n", argv[0]);
		exit(1);
	}
	// Using 1<<24 takes about 50 seconds to find gold.
	// int norder = 26;
	int norder = atoi(argv[1]);
	if (24 < norder)
		printf("Caution: large orders 24 < %d take a long time\n", norder);
	int nmax = (1<<norder) + 1;

	// Depth is how large any given sequence value can go.
	// Length is how long a sequence is.
	// Need to go to high depth to avoid big gap at golden mean.
	// int maxdepth = 16;
	// int maxlen = 9;
	int maxdepth = atoi(argv[2]);
	int maxlen = atoi(argv[3]);

	setup_gold(nmax);
	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	// iterate_fbaire(cfrac, 1, 10, 10, nmax);

	printf("#\n# Max order of polynomials = %d num=2^order = %d\n", norder, nmax);
	printf("#\n# Iterate to maxdepth=%d maxlen=%d\n#\n", maxdepth, maxlen);
	fflush (stdout);
	iterate_fbaire(cfrac, 1, maxdepth, maxlen, nmax);
#endif

// #define BINCOUNT_INDEX
#ifdef BINCOUNT_INDEX
	#define NBINS 337
	int gcount[NBINS];
	int fcount[NBINS];
	for (int i=0; i<NBINS; i++)
	{
		gcount[i] = 0;
		fcount[i] = 0;
	}

	int cfrac[SZ];
	long totg = 0;
	long totf = 0;

	// Max order 26 takes 6 seconds
	// 28 takes 20 seconds
	int maxord = 46;
	for (int ord=1; ord < maxord; ord++)
	{
		long gprev = totg;
		long nstart = 1UL << (ord-1);
		long nend = 2*nstart;
		double deno = (double) nstart;
		for (long n=nstart; n<nend; n ++)
		{
			int len = index_to_fbaire(cfrac, n);

			// Bin-count indexes
			double ex = ((double) (n-nstart)) / deno;
			ex *= NBINS;
			int bin = floor(ex + 0.5);

			// Reject the bad ones.
			if (len < 0)
			{
				fcount[bin]++;
				totf ++;
			}
			else
			{
				gcount[bin]++;
				totg ++;
printf("duude ord=%d n=%ld totg=%ld delta=%ld\n", ord+1, n, totg,
totg-gprev);
			}
		}
		printf("%d	%ld	%ld	%ld\n", ord+1, totg, totg-gprev, totf);
		fflush (stdout);
	}

#if 0
	printf("#\n# Total bincount = %d %d\n#\n", totg, totf);
	for (int i=0; i<NBINS; i++)
	{
		double ex = (((double) i) + 0.5) / ((double) NBINS);
		double pg = NBINS * ((double) gcount[i]) / ((double) totg);
		double pf = NBINS * ((double) fcount[i]) / ((double) totf);
		printf("%d	%g	%g	%g\n", i, ex, pg, pf);
	}
#endif
#endif

}
