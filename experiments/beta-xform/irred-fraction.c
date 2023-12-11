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


/* Return the beta value corresponding to the n'th golden polynomial.
 * It is be constructed from the bit string of (2n+1). Construction
 * is the mid-point construction: repeated iteration of the midpoint
 * 1/2 with this beta will (re-)generate the same bitstring, until
 * returning to the midpoint. The bits are just whether the orbit went
 * left or right of midpoint. The length of the orbit will be log_2(2n+1).
 */
double gold(unsigned long n, double x)
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
	double fmid = gold(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/** Helper array, needed for finding gold midpoints. */
double* zero = NULL;
void malloc_gold(int nmax)
{
	zero = (double*) malloc((nmax+1)*sizeof(double));
	for (int i=0; i<=nmax; i++) zero[i] = -1.0;

	// Allow valid memref to zero[-1] denoting beta=1.0
	zero++;
	zero[-1] = 1.0;
	zero[0] = 2.0;
}

/** Fill up array of zero candidates. Optionally needed. */
void fill_gold(long n)
{
	// Go top down, assuming lower ranks already filled.
	for (int i=n; 0<i; i--)
	{
		if (zero[i] > -0.5) break;
		zero[i] = find_zero(i, 1.0, 2.0);
	}
}

// Return true if the polynomial root is properly bracketed for
// the index specifying that polynomial.
bool zero_is_bracketed(int n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15
	// printf("---------\ncheck bracketing for gold=%20.16g at n=%d\n", gold, n);
	bool ork = true;
	int nhl = n;
	int nh = nhl >> 1;
	while (nh)
	{
		//printf("walk to n=%d nhl=%d nh=%d znh=%g go=%g comp=%d bad=%d\n",
		//       n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}
	// printf("Bracket says ork=%d\n", ork);

	return ork;
}

/**
 * Return the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not correspond to a golden polynomial,
 * return zero.
 */
double find_gold(long n)
{
	if (-1 == n) return 1.0;

	fill_gold(n);
	double gold = zero[n];

	if (zero_is_bracketed(n, gold)) return gold;
	return 0.0;
}

// =================================================================

// Perform standard iteration
double tee(double beta, double x)
{
	if (x<=0.5) return beta*x;
	return beta*(x-0.5);
}

// Compute the order of the iteration of the midpoint
// This is just sanity check; beta should "already be correct".
int iteration_order(double beta)
{
#define DELTA 1e-15
#define EPSI  1e-13
	int count = 0;
	double midpoint = 0.5-DELTA;
	do
	{
		count ++;
		midpoint = tee(beta, midpoint);
		// printf("%d  %g\n", count, midpoint);
	}
	while (fabs(midpoint-0.5)>EPSI && count < 100);
	return count;
}

// =================================================================

/** Return order of beta polynomial for index n.
 * This is just the length of binary bitstring for 2n+1.
 */
int order(unsigned long n)
{
	int len=0;
	unsigned long bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/**
 * Print index as golden bitsequence.
 * Actually print order, then index, then bitseq
 */
void prt_bitstr(unsigned long n, const char* pfx, const char* sfx)
{
	int ord = order(n);
	printf("%s o(%d) %ld={", pfx, ord, n);
	unsigned long bitstr = 2*n+1;
	for (int i=0; i<ord; i++)
	{
		printf("%ld", 0x1 & bitstr>>(ord-i-1));
	}
	printf("} %s", sfx);
}

// =================================================================

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

// =================================================================

/*
 * Given a finite-length Baire representation, return the
 * corresponding index number.
 */
long index_from_fbaire(int cfrac[], int len)
{
	// zero length corresponds to beta=1 which has index infinity
	// Which we report as -1;
	if (0 == len) return -1;

	if (1 == len)
	{
		// Sequence of [-1] is index zero which is beta=2
		if (-1 == cfrac[0]) return 0;
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

void print_seq(int cfrac[], int len, char* head, char* tail)
{
	printf("%s [", head);
	for (int i=0; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

// Generate finite-Baire sequence labels from a polynomial index.
// Given a polynomial index, set "cfrac" to the matching sequence.
// Return the length of the cfrac sequence.
// This is the inverse of what index_from_fbaire() does.
int index_to_fbaire(int cfrac[], unsigned long pidx)
{
	// #define DBG(X) printf X
	#define DBG(X)
	DBG(("enter index_to_fbaire pidx=%ld\n", pidx));

	// Count leading powers of two.
	int msum = 0;
	while (0 == pidx %2) { msum ++; pidx /=2; }

	// If it was a pure power, then terminate recursion.
	// This can only correspond to a sequence of length one,
	// the sequence labelled as [msum]. So set the digit
	// and return;
	if (1 == pidx)
	{
		DBG(("the end msum=%d\n", msum));
		cfrac[0] = msum;
		cfrac[1] = -777; // Terminator poison.
		return 1;
	}

	// If we are here, pidx is odd; truncate and try again.
	pidx = (pidx-1)/2;
	DBG(("red after odd to pidx=%ld msum=%d\n", pidx, msum));

	// Recurse.
	int len = index_to_fbaire(cfrac, pidx);
	DBG(("post recursion for pidx=%ld msum=%d new len=%d\n", pidx, msum, len));

	// Pass rejection slips back up the chain.
	if (len < 0) return len;

	// print_seq(cfrac, len, "to fb ", "\n");

	// What's left over is m_k
	int ord = order(pidx) - 2;
	DBG(("pre ord for pidx=%ld msum=%d ord=%d\n", pidx, msum, ord));
	msum -= ord;
	if (msum < 0) msum = 0;
	cfrac[len] = msum;
	cfrac[len+1] = -555; // Poison end-of-string marker

	// Return the length of the beast.
	return len+1;
}

// =================================================================
#define SZ 60

// Given polynomial index n, return the index defining the right
// side of the bracket containing this index.
long get_bracket_right(long n)
{
	int cfrac[SZ];
	int len = index_to_fbaire(cfrac, n);
	// printf("Enter get_right, n=%ld len=%d ", n, len);
	// print_seq(cfrac, len, "seq", "\n");

	int dig = len-1;

	// If the last digit is zero, truncate. Keep on truncating,
	// as long as we find zeros.
	if (0 == cfrac[dig])
	{
		do {
			len = dig;
			cfrac[dig] = -444;
			dig --;
		}
		while (0 <= dig && 0 == cfrac[dig]);
	}

	// Get the last digit. If it's positive, demote.
	if (0 <= dig && 0 < cfrac[dig])
		cfrac[dig]--;

	// Special case to convert empty sequence to [-1] denoting beta=2
	if (0 == len)
	{
		cfrac[0] = -1;
		len = 1;
	}

	// print_seq(cfrac, len, "right seq", "\n");
	long nright = index_from_fbaire(cfrac, len);
	// printf("right idx=%ld\n", nright);
	return nright;
}

// Given polynomial index n, return the index defining the left
// side of the bracket containing this index.
long get_bracket_left(long n)
{
	int cfrac[SZ];
	int len = index_to_fbaire(cfrac, n);
	// printf("Enter get_left, n=%ld len=%d ", n, len);
	// print_seq(cfrac, len, "seq", "\n");

	// Truncate the last digit.
	int dig = len-1;
	cfrac[dig] = -333;
	len = dig;

	long nleft = index_from_fbaire(cfrac, len);
	// printf("left idx=%ld\n", nleft);
	return nleft;
}

/*
 * Validate bracketing for betas and for the finite-Baire sequences.
 * Similar to validate_fbaire() below, but takes a more principled
 * approach.
 */
void validate_bracket(long n)
{
	double gold = find_gold(n);
	// If its not a valid index, do nothing.
	if (gold < 1.0) return;

	// printf("Validate bracket for %ld\n", n);
	// Verify gold bracketing

	// Verify the left bracket by ripping out powers of two until
	// an odd number is reached.
	long nleft = get_bracket_left(n);
	double gleft = find_gold(nleft);
	if (gleft < 0.5) printf("Error: no such left index %ld\n", n);
	if (gleft >= gold)
		printf("Error: bad left bracket at %ld: nr=%ld gold=%g gleft=%g\n",
			n, nleft, gold, gleft);

	// Verify right bracketing by knocking off only one power of two.
	long nright = get_bracket_right(n);
	double gright = find_gold(nright);
	if (gright < 0.5) printf("Error: no such right index %ld\n", (n-1)/2);
	if (gright <= gold)
		printf("Error: bad right bracket at %ld: nr=%ld gold=%g gright=%g\n",
			n, nright, gold, gright);
}

// =================================================================
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

// =================================================================

// #define VERFIY_FBAIRE
#ifdef VERFIY_FBAIRE

	// Verify reverse listing.
	// Takes 30 cpu-seconds to get to 1<<26
	int maxord = 26;
maxord=8;
	int toterr = 0;
	for (int ord=2; ord < maxord; ord++)
	{
		int totgood = 0;
		long nstart = 1UL << (ord-2);
		long nend = 2*nstart;
		for (long n=nstart; n<nend; n ++)
		{
			int cfrac[SZ];
			for (int i=0; i<SZ; i++) cfrac[i] = -666;
			int len = index_to_fbaire(cfrac, n);
			if (len < 0)
			{
#define PRINT_SEQS
#ifdef PRINT_SEQS
				printf(">>>>> %ld rejected\n", n);
#endif
				continue;
			}
			long seqno = index_from_fbaire(cfrac, len);
			if (n != seqno)
			{
				printf("Sequence numbering fail!! in=%ld out=%ld", n, seqno);
				toterr ++;
			}
#ifdef PRINT_SEQS
			printf(">>>>> %ld len=%d ", n, len);
			print_seq(cfrac, len, "sequence ", "\n");
#endif
			totgood ++;
		}
		printf("Observed %d betas at order %d\n", totgood, ord);
	}
	long nend = 1UL << maxord;
	printf("Verified up to order=%d n=%ld errors=%d\n", maxord, nend, toterr);
#endif

// #define REVERSE_VERIFY
#ifdef REVERSE_VERIFY
	// Verify that the reverse-indexing works.
	int norder = 20;
	int nmax = (1<<norder) + 1;

	malloc_gold(nmax);

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

	malloc_gold(nmax);
	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	// Iterating to length 10, depth 10 takes more than an hour,
	// mostly due to large numbers of overflow failures.
	// iterate_fbaire(cfrac, 1, 3, 8, nmax);
	iterate_fbaire(cfrac, 1, 3, 4, nmax);
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

	malloc_gold(nmax);
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

// =================================================================

int main(int argc, char* argv[])
{
#if 0
	int cfrac[SZ];
	for (int i=0; i<SZ; i++) cfrac[i] = -666;
	int nind = 53;
	// nind = 27;
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

// #define PRINT_INDEX
#ifdef PRINT_INDEX
	// Do nothing except print beta and indexes
	int nmax = 64;
	malloc_gold(nmax);
	for (int n=1; n<nmax; n ++)
	{
		double beta = find_gold(n);
		if (0.5 < beta)
		{
			printf("Good %d %g\n", n, beta);
		}
		else
		{
		}
	}
#endif

	int cfrac[SZ];
	malloc_gold(1024);
#if 0
	long ni = 53;
	double beta = find_gold(ni);
	int len = index_to_fbaire(cfrac, ni);
	print_seq(cfrac, len, "dbg", "\n");

	int order = iteration_order(beta);
	printf("got %d\n", order);
	prt_bitstr(ni, "bits", "\n");
#endif

	for (long n=1; n<16; n++)
	{
		double beta = find_gold(n);
		if (beta < 0.5) continue;

		validate_bracket(n);
		long nleft = get_bracket_left(n);
		long nright = get_bracket_right(n);
		printf("bracket (%ld |=> %ld <=| %ld) ", nleft, n, nright);
		printf("gold=%g ", beta);
		prt_bitstr(n, "bits", "");

		for (int i=0; i<SZ; i++) cfrac[i] = -666;
		int len = index_to_fbaire(cfrac, n);
		if (len < 0)
			printf("\nError: missing representation for n=%ld\n", n);
		print_seq(cfrac, len, "", "\n");
		printf("-------\n");
	}

}
