/*
 * irred-fraction.c
 *
 * Find integer sequences for the golden polynomials.
 * Relate it to the continued-fraction representation.
 *
 * February 2018, October 2020, December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "irred-gold.c"

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
	if (0 == len) printf("INF");
	if (0 < len) printf("%d", cfrac[0]);
	for (int i=1; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

/*
 * Given a finite-length Baire representation, return the
 * corresponding index number.
 */
long index_from_fbaire(int cfrac[], int len)
{
	#define DST(X)
	// #define DST(X) X
	DST(print_seq(cfrac, len, "enter index_from_fbaire", "\n"));

	if (0 > len)
	{
		printf("Mega failure len=%d\n", len);
		return -666;
		// *((int*) 0x42) = 33;
	}

	// zero length corresponds to beta=1 which has index infinity
	// Which we report as -1;
	if (0 == len) return -1;

	if (1 == len)
	{
		// Sequence of [-1] is index zero which is beta=2
		if (-1 == cfrac[0]) return 0;
		return 1UL << cfrac[0];
	}

	long bracket = index_from_fbaire(cfrac, len-1);
	if (0 > bracket) return bracket; // report overflow

	long leader = find_leader(bracket);
	DST(printf("leader is %ld\n", leader));
	if (0 > leader) return leader; // report overflow

#if OLD_BRUTE_FORCE
	long leader = 2UL * bracket + 1UL;
	DST(printf("leader is %ld\n", leader));

	// Cannot issue a sequence number that does not correspond to
	// a valid golden polynomial. I can't guess what's valid or not,
	// so go to the polynomial itself for ground truth.
	int cnt = 0;
	while (false == is_valid_index(leader) && cnt < 60)
	{
		leader *= 2UL;
		cnt ++;
	}
	if (60 < cnt + cfrac[len-1])
	{
		// printf("Error: Overflow during sequence decode\n");
		return -444;
	}
#endif

	// Trailing digit encodes index-doubling
	leader *= 1UL << cfrac[len-1];

	DST(printf("exit index_from_baire leader is %ld after shift=%d\n", leader, shift));
	return leader;
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

	// Index of -1 maps to beta=1 and the empty sequence
	if (-1 == pindex)
		return 0;

	// Index of zero maps to beta=2 so [-1] the infinite sequence of zeros.
	if (0 == pindex)
	{
		cfrac[0] = -1;
		return 1;
	}

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
	if (0 > base) return base; // overflow condition

	pidx = pindex / base;
	if (0 == pidx) {
		printf("FATAL algo error!\n");
		pidx = 1; // terrible algorithm fail
	}

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
	// print_seq(cfrac, len, "seq ", "\n");

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

	// print_seq(cfrac, len, "computed right seq ", "\n");
	long nright = index_from_fbaire(cfrac, len);
	// printf("computed right idx=%ld\n", nright);
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

	// length of zero corresponds to beta=1 which is mapped to index -1
	if (0 == len) return -1;

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

	if (-1 > n)
		{ printf("Error: overflow index %ld\n", n); return false; }

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
	if (gleft >= gold && -1 != n)
		{ printf("Error: bad left bracket at %ld: nleft=%ld gold=%g gleft=%g\n",
			n, nleft, gold, gleft); ok = false; }

#if 1
	// Verify right bracketing by knocking off only one power of two.
	long nright = get_bracket_right(n);
	double gright = find_gold(nright);
	if (gright < 0.5)
	{
		printf("Error: no such right index %ld for %ld\n", nright, n);
		int cfrac[SZ];
		int len = index_to_fbaire(cfrac, n);
		print_seq(cfrac, len, "Index has seq ", "");
		len = index_to_fbaire(cfrac, nright);
		print_seq(cfrac, len, " <=| bad right seq ", "\n");
		ok = false;
	}
	if (gright <= gold)
		{ printf("Error: bad right bracket at %ld: nright=%ld gold=%g gright=%g\n",
			n, nright, gold, gright); ok = false; }
#endif

	// Validate conversion to and from Baire.
	int cfrac[SZ];
	int len = index_to_fbaire(cfrac, n);
	long seqno = index_from_fbaire(cfrac, len);
	if (n != seqno)
	{
		printf("Sequence numbering fail!! in=%ld out=%ld ", n, seqno);
		print_seq(cfrac, len, "seq ", "\n");
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
	printf("%ld	%g	%g	%g %ld	%g	%ld	%g # ",
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
	// To avoid having everything bunch up near beta=2,
	// also decrease the depth as length increases.
	if (len < maxlength)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len] = 0;
		recurse_fbaire(bfrac, len+1, maxdepth/2, maxlength, maxn, do_print);
	}

	long idx = index_from_fbaire(cfrac, len);

	// Don't bother with validation if out of bounds
	if (idx < maxn && -1L < idx)
	{
		validate_bracket(idx);

		// Print equivalent continued fraction, for the odometer graph
		if (do_print) print_odo_graph(cfrac, len);
	}

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
	long nmax = (1UL<<norder);

	malloc_gold(nmax+1);
	fill_gold(nmax);

	int cfrac[SZ];
	cfrac[0] = 0;
	recurse_fbaire(cfrac, 1, depth, length, nmax, do_print);
}

// =================================================================

void bincount_index()
{
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
}

void print_debug_info(long seqno)
{
	malloc_gold(seqno+1);

	int cfrac[SZ];
	int slen = index_to_fbaire(cfrac, seqno);
	printf("Index: %ld len=%d", seqno, slen);
	print_seq(cfrac, slen, " sequence ", "\n");

	double beta = find_gold(seqno);

	if (0.5 < beta) printf("beta = %g\n", beta);

	if (beta < 0.5)
	{
		beta = find_poly_zero(seqno);
		long factor = theta_factor(seqno, beta);
		double feta = find_poly_zero(factor);
		slen = index_to_fbaire(cfrac, factor);
		printf("Invalid index; factors to %g = %ld = ", feta, factor);
		print_seq(cfrac, slen, "", "\n");
	}

#if 1
	validate_bracket(seqno);

	long nleft = get_bracket_left(seqno);
	double gleft = find_gold(nleft);
	slen = index_to_fbaire(cfrac, nleft);
	printf("Left limit: %ld = %g = ", nleft, gleft);
	print_seq(cfrac, slen, "", "\n");

	long nright = get_bracket_right(seqno);
	double gright = find_gold(nright);
	slen = index_to_fbaire(cfrac, nright);
	printf("Right limit: %ld = %g = ", nright, gright);
	print_seq(cfrac, slen, "", "\n");
#endif
}

// =================================================================

int main(int argc, char* argv[])
{
	// bincount_index();

// #define SEQUENCE_EXPLORER
#ifdef SEQUENCE_EXPLORER
	// Obtain one sequence from command line. Print it's index.
	if (1 == argc) {
		fprintf(stderr, "Usage: %s <sequence>\n", argv[0]);
		exit(1);
	}
	int len = argc-1;
	int cfrac[SZ];
	for (int i=0; i<len; i++) cfrac[i] = atoi(argv[i+1]);
	long seqno = index_from_fbaire(cfrac, len);

	printf("Index: %ld len=%d", seqno, len);
	print_seq(cfrac, len, " provided ", "\n");

	print_debug_info(seqno);
#endif

// #define INDEX_EXPLORER
#ifdef INDEX_EXPLORER
	// Obtain one index from command line. Print debug info for it.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}
	long nidx = atol(argv[1]);
	print_debug_info(nidx);
#endif

// #define PRINT_INDEX
#ifdef PRINT_INDEX
	// Do nothing except print beta and indexes
	// Obtain max to print from command line.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <max-index>\n", argv[0]);
		exit(1);
	}
	long nmax = atol(argv[1]);
	malloc_gold(nmax);
	fill_gold(nmax);
	for (long n=1; n<nmax; n ++)
	{
		double beta = find_gold(n);
		if (0.5 < beta)
		{
			printf("Good %ld %g\n", n, beta);
		}
		else
		{
		}
	}
#endif

// #define VALIDATE_INDEX
#ifdef VALIDATE_INDEX
	// Validate indexes in sequential order.
	// Obtain max index from command line.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <max-order>\n", argv[0]);
		exit(1);
	}
	int nord = atoi(argv[1]);
	long nmax = 1UL << nord;
	malloc_gold(nmax);

	int cfrac[SZ];
	for (long n=1; n<nmax; n++)
	{
		double beta = find_gold(n);
		if (beta < 0.5) continue;

		validate_bracket(n);
		for (int i=0; i<SZ; i++) cfrac[i] = -666;
		int len = index_to_fbaire(cfrac, n);
		if (len < 0)
			printf("\nError: missing representation for n=%ld\n", n);

		// if (4 == len && 0 == cfrac[3] && 0 != cfrac[0])
		if (1)
		{
			long nleft = get_bracket_left(n);
			long nright = get_bracket_right(n);
			printf("bracket (%ld |=> %ld <=| %ld) ", nleft, n, nright);
			printf("=%g=", beta);
			prt_bitstr(n, " ", "");
			print_seq(cfrac, len, " ", "\n");
			// printf("-------\n");
		}
	}
#endif

// #define RECURSIVE_CHECK
#ifdef RECURSIVE_CHECK
	// Run validation on the recursively-generated sequences.
	// Same as the odometer graph below, but does not print data.
	if (4 != argc) {
		fprintf(stderr, "Usage: %s <order> <maxdepth> <maxlen>\n", argv[0]);
		exit(1);
	}

	int norder = atoi(argv[1]);
	int maxdepth = atoi(argv[2]);
	int maxlen = atoi(argv[3]);

	// Using order==24 takes about 50 seconds to find gold.
	if (24 < norder)
		printf("Caution: large orders 24 < %d take a long setup time\n", norder);
	if (60 < norder) exit(-1);

	double space = pow (maxdepth, maxlen);
	long niter = (long) floor(space);
	if ((1UL<<norder) < niter) niter = 1UL<<norder;
	printf("Predict %ld iterations\n", niter);
	if (1UL<<44 < niter) exit(-1);

	generate_fbaire(norder, maxdepth, maxlen, false);
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

	// Depth is how large any given sequence value can go.
	// Length is how long a sequence is.
	// Need to go to high depth to avoid big gap at golden mean.
	int maxdepth = atoi(argv[2]);
	int maxlen = atoi(argv[3]);

	// Using 1<<24 takes about 50 seconds to find gold.
	int norder = atoi(argv[1]);
	if (24 < norder)
		printf("Caution: large orders 24 < %d take a long setup time\n", norder);
	if (60 < norder) exit(-1);

	double space = pow (maxdepth, maxlen);
	long niter = (long) floor(space);
	if ((1UL<<norder) < niter) niter = 1UL<<norder;
	printf("Predict %ld iterations\n", niter);
	if (1UL<<44 < niter) exit(-1);

	long nmax = (1UL<<norder) + 1;

	printf("#\n# Max order of polynomials = %d num=2^order = %ld\n", norder, nmax);
	printf("#\n# Iterate to maxdepth=%d maxlen=%d\n", maxdepth, maxlen);
	printf("# Predict %ld iterations\n", niter);
	printf("#\n");
	fflush (stdout);
	generate_fbaire(norder, maxdepth, maxlen, true);
#endif
}
