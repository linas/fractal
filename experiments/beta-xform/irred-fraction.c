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
long zero_bracket_factor(long n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15
	// printf("---------\ncheck bracketing for gold=%20.16g at n=%d\n", gold, n);
	bool ork = true;
	long nhl = n;
	long nh = nhl >> 1;
	while (nh)
	{
		//printf("walk to n=%ld nhl=%ld nh=%ld znh=%g go=%g comp=%d bad=%d\n",
		//       n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}
	// printf("Bracket says ork=%d\n", ork);

	if (ork) return -1;
	return nh;
}

/**
 * Return the single, unique real zero of the n'th bitstring polynomial.
 * Recall most values of `n` do NOT correspond to a golden polynomial.
 */
double find_poly_zero(long n)
{
	fill_gold(n);
	return zero[n];
}

/**
 * Return the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not correspond to a golden polynomial,
 * return zero.
 */
double find_gold(long n)
{
	if (-1 == n) return 1.0;
	double gold = find_poly_zero(n);
	if (-1 == zero_bracket_factor(n, gold)) return gold;
	return 0.0;
}

// =================================================================

// Perform standard iteration
double tee(double beta, double x)
{
	if (x<=0.5) return beta*x;
	return beta*(x-0.5);
}

// Compute the period (order) of the iteration of the midpoint
// This is just sanity check; beta should "already be correct".
int iteration_period(double beta)
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
	printf("%s|%d| %ld={", pfx, ord, n);
	unsigned long bitstr = 2*n+1;
	for (int i=0; i<ord; i++)
	{
		printf("%ld", 0x1 & bitstr>>(ord-i-1));
	}
	printf("}%s", sfx);
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

void print_seq(int cfrac[], int len, char* head, char* tail)
{
	printf("%s[", head);
	if (0 <= len) printf("%d", cfrac[0]);
	for (int i=1; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

/*
 * Given a finite-length Baire representation, return the
 * corresponding index number.
 */
long index_from_fbaire(int cfrac[], int len)
{
	// print_seq(cfrac, len, "enter index_from_fbaire", "\n");

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
	// printf("leader is %ld\n", follower);

	// Trailing digit encodes index-doubling
	int shift = cfrac[len-1];

	// Leading digit provides unique coding for the 2[]+1 operation
	if (1<len) shift += cfrac[0];

#if 0
	int bump = len-3;
	if (bump < 0) bump = 0;
	for (int j=0; j<= bump; j++)
		shift += cfrac[j];
printf("shift=%d\n", shift);
#endif

	// More careful overflow check
	int nbits = 0;
	long fo = follower;
	while (fo >>= 1) ++nbits;
	if (60 < nbits+shift) return -1;

	follower *= 1UL << shift;

	return follower;
}

// Generate finite-Baire sequence labels from a polynomial index.
// Given a polynomial index, set "cfrac" to the matching sequence.
// Return the length of the cfrac sequence.
// This is the inverse of what index_from_fbaire() does.
int index_to_fbaire(int cfrac[], unsigned long pindex)
{
	// #define DBG(X) printf X
	#define DBG(X)
	DBG(("enter index_to_fbaire pidx=%ld\n", pindex));

	// Count leading powers of two.
	long pidx = pindex;
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

	// After recursion, we always append a new digit.
	// Start by appending a zero, then get what the resulting
	// sequence encodes. The ratio of that, to the pindex should
	// be an exact power of two. Get that, and that will be our
	// trailing digit.  The need for this is to distinguish
	// index 10 which is a leader, but is a muliple of (invalid) index 5
	// and 14, which is a follower of 7.
	cfrac[len] = 0;
	cfrac[len+1] = -555; // Poison end-of-string marker
	// print_seq(cfrac, len+1, "baseline ", "\n");
	long base = index_from_fbaire(cfrac, len+1);
	pidx = pindex / base;
if (0 == pidx) pidx = 1; // fail

	msum = 0;
	while (0 == pidx %2) { msum ++; pidx /=2; }
	cfrac[len] = msum;
	DBG(("post-base for pidx=%ld base=%ld msum=%d \n", pindex, base, msum));

	DBG(("-----\n"));
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
 */
bool validate_bracket(long n)
{
	bool ok = true;
	double gold = find_gold(n);

	if (gold < 1.0)
		{ printf("Error: invalid index %ld\n", n); ok = false; }

	// printf("Validate bracket for %ld\n", n);
	// Verify gold bracketing

	// Verify the left bracket by ripping out powers of two until
	// an odd number is reached.
	long nleft = get_bracket_left(n);
	double gleft = find_gold(nleft);
	if (gleft < 0.5)
		{ printf("Error: no such left index %ld for %ld\n", nleft, n); ok = false; }
	if (gleft >= gold)
		{ printf("Error: bad left bracket at %ld: nr=%ld gold=%g gleft=%g\n",
			n, nleft, gold, gleft); ok = false; }

	// Verify right bracketing by knocking off only one power of two.
	long nright = get_bracket_right(n);
	double gright = find_gold(nright);
	if (gright < 0.5)
		{ printf("Error: no such right index %ld for %ld\n", nright, n); ok = false; }
	if (gright <= gold)
		{ printf("Error: bad right bracket at %ld: nr=%ld gold=%g gright=%g\n",
			n, nright, gold, gright); ok = false; }

	// Validate conversion to and from Baire.
	int cfrac[SZ];
	int len = index_to_fbaire(cfrac, n);
	long seqno = index_from_fbaire(cfrac, len);
	if (n != seqno)
	{
		printf("Sequence numbering fail!! in=%ld out=%ld ", n, seqno);
		print_seq(cfrac, len, "seq", "\n");
		ok = false;
	}
	if (!ok) printf("-------\n");

	return ok;
}

// =================================================================
/*
 * Given label sequence, print equivalent continued-fraction value.
 * Used for producing odometer graph.
 */
void print_odo_graph(int cfrac[], int len)
{
	long idx = index_from_fbaire(cfrac, len);
	double gold = find_gold(idx);

	long nleft = get_bracket_left(idx);
	double gleft = find_gold(nleft);

	long nright = get_bracket_right(idx);
	double gright = find_gold(nright);

	double ex = cf_to_real(cfrac, len);
	double rex = reverse_cf_to_real(cfrac, len);
	printf("%ld	%g	%g	%g %ld	%g	%ld	%g #",
		idx, gold, ex, rex, nleft, gleft, nright, gright);

	int dfrac[SZ];
	int dlen = index_to_fbaire(dfrac, nleft);
	print_seq(dfrac, dlen, "", " |=> ");
	print_seq(cfrac, len, "", "");

	dlen = index_to_fbaire(dfrac, nright);
	print_seq(dfrac, dlen, " <=| ", "\n");
	fflush(stdout);
}

/*
 * Generate correctly-ordered finite-baire sequences. The ordering
 * is such that the corresponding golden number is strictly decreasing
 * for the generated cfrac sequences.
 *
 * maxdepth == max number of doubling steps.
 * maxlength == max length of sequence.
 * maxn == cutoff for highest known n
 * do_print == boolean, print the sequences
 */
void recurse_fbaire(int cfrac[], int len,
                    int maxdepth, int maxlength, long maxn,
                    bool do_print)
{
	// Iterate to max length, first.
	if (len < maxlength)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len] = 0;
		recurse_fbaire(bfrac, len+1, maxdepth, maxlength, maxn, do_print);
	}

	long idx = index_from_fbaire(cfrac, len);
	if (idx >= maxn) return; // why ??
	validate_bracket(idx);

	// Print equivalent continued fraction, for the odometer graph
	if (do_print) print_odo_graph(cfrac, len);

	// Iterate depthwise second.
	if (cfrac[len-1] < maxdepth)
	{
		int afrac[SZ];
		for (int i=0; i<len; i++) afrac[i] = cfrac[i];
		afrac[len-1] ++;
		recurse_fbaire(afrac, len, maxdepth, maxlength, maxn, do_print);
	}
}

/*
 * Generate correctly-ordered finite-baire sequences. The ordering
 * is such that the corresponding golden number is strictly decreasing
 * for the generated cfrac sequences.
 *
 * norder == max polynomial order to go to (order == length of bitstring/orbit)
 * maxdepth == max number of index doubling steps
 * maxlength == max length of sequence.
 * do_print == boolean, print the sequences
 */
void generate_fbaire(int norder, int depth, int length, bool do_print)
{
	int nmax = (1<<norder) + 1;

	malloc_gold(nmax);
	fill_gold(nmax);

	int cfrac[SZ];
	cfrac[0] = 0;
	recurse_fbaire(cfrac, 1, depth, length, nmax, do_print);
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
	print_seq(cfrac, len, " provided ", "\n");

	int slen = index_to_fbaire(cfrac, seqno);
	printf("Index: %d len=%d", seqno, slen);
	print_seq(cfrac, slen, " reconstruct ", "\n");

	malloc_gold(seqno+1);
	double beta = find_gold(seqno);

	if (0.5 < beta) printf("beta = %g\n", beta);

	if (beta < 0.5)
	{
		beta = find_poly_zero(seqno);
		long factor = zero_bracket_factor(seqno, beta);
		double feta = find_poly_zero(factor);
		slen = index_to_fbaire(cfrac, factor);
		printf("Invalid index; factors to %g = %ld = ", feta, factor);
		print_seq(cfrac, slen, "", "\n");
	}
#endif

// #define PRINT_INDEX
#ifdef PRINT_INDEX
	// Do nothing except print beta and indexes
	int nmax = 64;
	malloc_gold(nmax);
	fill_gold(nmax);
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

#ifdef NEW_STUFF
	int cfrac[SZ];
#if 0
	malloc_gold(100);
	long ni = 53;
	double beta = find_gold(ni);
	int len = index_to_fbaire(cfrac, ni);
	print_seq(cfrac, len, "dbg", "\n");

	int order = iteration_period(beta);
	printf("got %d ", order);
	prt_bitstr(ni, "bits", "\n");
	exit(0);
#endif

	int nord = 1<<16;
	malloc_gold(nord);
	for (long n=1; n<nord; n++)
	{
		double beta = find_gold(n);
		if (beta < 0.5) continue;

		validate_bracket(n);
		for (int i=0; i<SZ; i++) cfrac[i] = -666;
		int len = index_to_fbaire(cfrac, n);
		if (len < 0)
			printf("\nError: missing representation for n=%ld\n", n);

		if (4 == len && 0 == cfrac[3] && 0 != cfrac[0])
		{
			long nleft = get_bracket_left(n);
			long nright = get_bracket_right(n);
			printf("bracket (%ld |=> %ld <=| %ld) ", nleft, n, nright);
			printf("=%g=", beta);
			prt_bitstr(n, "", "");
			print_seq(cfrac, len, "", "\n");
			printf("-------\n");
		}
	}
#endif

// #define SANITY_CHECK
#ifdef SANITY_CHECK
	// Run validation on the recursively-generated sequences.
	// Same as the odometer graph below, but does not print data.
	if (4 != argc) {
		fprintf(stderr, "Usage: %s <order> <maxdepth> <maxlen>\n", argv[0]);
		exit(1);
	}

	int norder = atoi(argv[1]);
	int maxdepth = atoi(argv[2]);
	int maxlength = atoi(argv[3]);

	// Using order==24 takes about 50 seconds to find gold.
	if (24 < norder)
		printf("Caution: large orders 24 < %d take a long time\n", norder);

	if (10 < maxdepth)
		printf("Caution: large depths 10 < %d take a long time\n", maxdepth);

	generate_fbaire(norder, maxdepth, maxlength, false);
#endif

#define ODOMETER_GRAPH
#ifdef ODOMETER_GRAPH
	// Generate expansions in sequential order, then print the equivalent
	// index, beta and continued-frac equivalent. Used to make the odometer
	// graph for the paper.
	if (4 != argc) {
		fprintf(stderr, "Usage: %s <order> <maxdepth> <maxlen>\n", argv[0]);
		exit(1);
	}

	// Using 1<<24 takes about 50 seconds to find gold.
	int norder = atoi(argv[1]);
	if (24 < norder)
		printf("Caution: large orders 24 < %d take a long time\n", norder);

	// Depth is how large any given sequence value can go.
	// Length is how long a sequence is.
	// Need to go to high depth to avoid big gap at golden mean.
	int maxdepth = atoi(argv[2]);
	int maxlen = atoi(argv[3]);

	if (10 < maxdepth)
		printf("Caution: large depths 10 < %d take a long time\n", maxdepth);

	int nmax = (1<<norder) + 1;

	printf("#\n# Max order of polynomials = %d num=2^order = %d\n", norder, nmax);
	printf("#\n# Iterate to maxdepth=%d maxlen=%d\n#\n", maxdepth, maxlen);
	fflush (stdout);
	generate_fbaire(norder, maxdepth, maxlen, true);
#endif

}
