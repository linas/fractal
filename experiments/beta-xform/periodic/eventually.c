/*
 * eventually.c
 *
 * Code for eventually-periodic sequences
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"
#include "selfie-rational.c"
#include "selfie-tree.c"

// ---------------------------------------------------------------------
/* Check for prefix+cycle combinations that are not minimal.
 * Returns true if in minimal form, else return false.
 *
 * Checks several conditions:
 * -- Does the prefx+cycle combo have the shortest possible prefix?
 * -- Is the cycle minimal
 */
bool is_prefix_ok(long pfx, long cyc, int cyclen)
{
	if (0 == pfx) return true;

	int pfxlen = bitlen(pfx);

	// If prefix identical to cycle... then prefix can be absorbed into
	// the cycle to get a shorter expression.
	if (pfxlen == cyclen && pfx == cyc) return false;

	// Check cyclic permutation = if last bit of prefix can be rotated
	// into the cyclic sequence, then prefix is not maximally short.
	if ((pfx & 1UL) == (cyc & 1UL)) return false;

	// Unexpected thing
	short cof[MAX_EVENT_COF];
	get_event_coeffs(cof, pfx, cyc, cyclen);
	int clen = 0;
	for (int i=0; -100 < cof[i]; i++) clen++;
	if (0 == cof[clen-1])
	{
		printf("Not minimal at (%ld, %ld/%d)\n", pfx, cyc, cyclen);
		return false;
	}

	// Make sure the cycle is minimal.
	// If cyclen is a multiple of a prime, check for equivalent shorter lengths.
	// Checking multiples of 2 and 3 are enough for our purposes.
	if (3< cyclen && 0 == cyclen %2)
	{
		long a = cyc & ((1UL<< (cyclen/2)) -1UL);
		long b = cyc >> (cyclen/2);
		if (a==b) return false;
	}

	if (5< cyclen && 0 == cyclen %3)
	{
		long mask = (1UL<< (cyclen/3)) - 1UL;
		long a = cyc & mask;
		long b = (cyc >> (cyclen/3)) & mask;
		long c = cyc >> (2*cyclen/3);
		if (a==b && a==c) return false;
	}

	return true;
}

/*
 * Rotate the cycle left, until the first bit of the cycle is equal to 'bit'.
 * 'bit' can be zero or one.
 * The prefix is extended as the rotations is done.
 * Return true if changed.
 */
bool rotate_until_bit(long* pfxp, long* cycp, int cyclen, int bit)
{
	// if (0 == cyclen) return false;
	// if (0 == *cycp) return false;
	long pfx = *pfxp;
	long cyc = *cycp;
	long hibit = 1UL << (cyclen-1);
	bit <<= (cyclen-1);
	if (bit == (cyc & hibit)) return false;

	while (bit != (cyc & hibit))
	{
		pfx <<= 1;
		if (0 == bit) pfx |= 1UL;
		cyc <<= 1;
		if (0 == bit) cyc |= 1UL;
	}

	*pfxp = pfx;
	*cycp = cyc;
	return true;
}

/* Return true if the non-periodic version is valid. */
bool is_valid_finite(long pfx, long cyc, int cyclen)
{
	rotate_until_bit(&pfx, &cyc, cyclen, 1);
	cyc |= 1UL;
	pfx <<= cyclen;
	pfx |= cyc;
	long idx = (pfx - 1UL) / 2;
	return valid_gold_index(idx);
}

#if LATER
/* Interpret the orbit as being a dyadic bitstring, and return
 * the corresponding double-precision value. Result will be between
 * zero and one.
 */
long event_to_bitseq(long pfx, long cyc, int cyclen)
{
	// Brute force, no finesse.
	int maxlen = 56;
	if (62-cyclen < maxlen) maxlen = 62 - cyclen;

	int len = bitlen(pfx);
	while (len < 62-cyclen)
	{
		pfx <<= cyclen;
		pfx |= cyc;
		len += cyclen;
	}

	// Use the same #define MAXBITS 60 as in orbit_to_bitseq
	if (MAXBITS < len) pfx >>= len-MAXBITS;
	if (MAXBITS > len) pfx <<= MAXBITS-len;

	return pfx;
}

double event_to_double(long pfx, long cyc, int cyclen)
{
	// Brute force, no finesse.
	int maxlen = 56;
	if (62-cyclen < maxlen) maxlen = 62 - cyclen;

	int len = bitlen(pfx);
	while (len < 62-cyclen)
	{
		pfx <<= cyclen;
		pfx |= cyc;
		len += cyclen;
	}
	unsigned long deno = 1UL << len;
	double dyd = ((double) pfx) / ((double) deno);
	dyd = 2.0 * (dyd - 0.5);
	return dyd;
}

/* Convert the finite verion of periodic orbit to double. */
double finite_to_double(long pfx, long cyc, int cyclen)
{
	rotate_until_bit(&pfx, &cyc, cyclen, 1);
	cyc |= 1UL;
	pfx <<= cyclen;
	pfx |= cyc;
	int len = bitlen(pfx);
	unsigned long deno = 1UL << len;
	double dyd = ((double) pfx) / ((double) deno);
	dyd = 2.0 * (dyd - 0.5);
	return dyd;
}
#endif

// cheap hack
void low_guess(double rat, int* p, int* q)
{
	for (int i=2; i<16386; i++)
	{
		if ((fabs(fmod(rat*i, 1.0)) < 1e-6) ||
		    (fabs(1.0-fmod(rat*i, 1.0)) < 1e-6))
		{
			*q = i;
			*p = floor(rat*i);
			return;
		}
	}
}

// ---------------------------------------------------------------------

void print_debug_info(unsigned long p, unsigned long q)
{
	unsigned long pfx;
	unsigned long cyc;
	int cyclen;
	get_event_cycle (p, q, &pfx, &cyc, &cyclen);
	printf("Prefix = %ld cycle=%ld len=%d\n", pfx, cyc, cyclen);

	bool goodpf = is_prefix_ok(pfx, cyc, cyclen);
	if (goodpf)
		printf("Prefix looks good\n");
	else
		printf("Error: prefix is reducible!!!\n");

	short cof[MAX_EVENT_COF];
	get_event_coeffs(cof, pfx, cyc, cyclen);
	print_event_coeffs(cof, "Coeffs are: ", "\n");

	double gold = find_event_zero(cof, 1.0, 2.0);
	printf("Root %20.16g\n", gold);

	bool pok = is_event_ok(p, q);
	if (pok)
		printf("Seems to be a valid polynomial\n");
	else
		printf("Error: Rejected polynomial\n");

	unsigned long dyad = beta_to_dyadic(gold);
	unsigned long drat = rational_to_dyadic(p+q, 2*q, WORDLEN);
	int badbit = compare_dyadics(dyad, drat);
	if (0 <= badbit && badbit < 40)
		printf("Error: bad orbit at %d\n", badbit);

	print_dyadic(dyad, 60, "midpoint orbi= ", "\n");
	print_dyadic(drat, 60, "rational dyad= ", "\n");
	if (-1 < badbit)
		printf("First miscompare at %d\n", badbit);
	else
		printf("Perfect compare!\n");

	unsigned long idx = beta_to_index(gold);
	printf("Index: %ld\n", idx);
	unsigned long moves = idx_to_moves(idx);
	print_moves(moves, "Tree moves: ", "\n");

	unsigned long pm, qm;
	moves_to_rational(moves, &pm, &qm);
	double mrat = ((double) pm) / ((double) qm);
	printf("Move rational =%ld/%ld  = %g\n", pm, qm, mrat);

	int pg, qg;
	low_guess(mrat, &pg, &qg);
	printf("Reduced guess = %d/%d\n", pg, qg);
}

int lesser(const void * px, const void * py)
{
	double x = *(double *) px;
	double y = *(double *) py;
	return x>y;
}

int main(int argc, char* argv[])
{
#define MANUAL_EXPLORE
#ifdef MANUAL_EXPLORE
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s <p> <q>\n", argv[0]);
		exit(1);
	}

	long p = atol(argv[1]);
	long q = atol(argv[2]);

	print_debug_info(p, q);
#endif

// #define BULK_LISTING
#ifdef BULK_LISTING
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <maxdeno>\n", argv[0]);
		exit(1);
	}

	long maxdeno = atol(argv[1]);

#define NROOTS 1000
	double roots[NROOTS];
	int totfnd = 0;
	for (long deno = 2; deno <= maxdeno; deno++)
	{
		for (long num = 1; num <deno; num++)
		{
			// Work only with reduced rationals
			if (1 < gcd(num, deno)) continue;

			// Don't look at those with zero prefix length
			// (These are the finite orbits)
			unsigned long pfx, cyc;
			int clen;
			get_event_cycle(num, deno, &pfx, &cyc, &clen);
			if (0 == pfx) continue;

// #define TEST_HYPOTHESIS
#ifdef TEST_HYPOTHESIS
			// Hypothesis testing. Set hypo see how it predicts.
			bool hypo = (pfx%2 == 1);  // False, see 7/12 and 13/20
			// bool hypo = (pfx != 2);    // true, never happens
			// int ell = bitlen(pfx);
			// bool hypo = (1==ell) || (pfx != (1UL<<(ell-1)));   // true.
			// bool hypo = (0 == cyc) || (pfx != 5);  // False; see 9/28

			bool evok = is_event_ok(num, deno);
			if (!hypo && evok)
			{
				printf(">>>> Reject valid  at  p/q = %ld/%ld\n", num, deno);
				print_debug_info(num, deno);
			}
			if (hypo && !evok)
			{
				// printf(">>>> Accept invalid  at  p/q = %ld/%ld\n", num, deno);
				// print_debug_info(num, deno);
			}
#endif

			if (false == is_event_ok(num, deno))
			{
				printf("fail at p/q = %ld/%ld\n", num, deno);
				// print_debug_info(num, deno);
				continue;
			}

			double gold = event_gold(num, deno);
			roots[totfnd] = gold;
			totfnd ++;
		}
		printf("------ done with %ld\n", deno);
	}
	printf("found %d\n", totfnd);

	// Verify that the roots really are distinct
	qsort(roots, totfnd, sizeof(double), lesser);
	for (int i=0; i<totfnd-1; i++)
	{
		// printf("%d %g\n", i, roots[i]);
		if (roots[i+1]-roots[i] < 1e-6)
			printf("Error: oops %d %g\n", i, roots[i]);
	}
#endif

// #define GRAPH_LISTING
#ifdef GRAPH_LISTING
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <maxdeno>\n", argv[0]);
		exit(1);
	}
	long maxdeno = atol(argv[1]);

	for (long deno = 2; deno <= maxdeno; deno++)
	{
		for (long num = 1; num <deno; num++)
		{
			// Work only with reduced rationals
			if (1 < gcd(num, deno)) continue;

			// Don't look at those with zero prefix length
			// (These are the finite orbits)
			unsigned long pfx, cyc;
			int clen;
			get_event_cycle(num, deno, &pfx, &cyc, &clen);
			if (0 == pfx) continue;

			bool isok = is_event_ok(num, deno);
			double gold = event_gold(num, deno);
			double ex = ((double) num) / ((double) deno);

			printf("%ld	%ld	%g	%d	%g\n", num, deno, ex, isok, gold);
		}
		printf("\n");
	}
#endif

// #define COMB
#ifdef COMB
	// This draws vertical tick-marks at the location of the valid periodic
	// sequences. This confirms the structure of the trimmed tree.
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <maxdeno>\n", argv[0]);
		exit(1);
	}
	long maxdeno = atol(argv[1]);

	for (long deno = 2; deno <= maxdeno; deno++)
	{
		for (long num = 1; num <deno; num++)
		{
			// Work only with reduced rationals
			if (1 < gcd(num, deno)) continue;

			// Don't look at those with zero prefix length
			// (These are the finite orbits)
			unsigned long pfx, cyc;
			int clen;
			get_event_cycle(num, deno, &pfx, &cyc, &clen);
			if (0 == pfx) continue;

			bool isok = is_event_ok(num, deno);
			double gold = event_gold(num, deno);
			double ex = ((double) num) / ((double) deno);

			// Draw tick-marks so we can see where they are.
			if (isok)
			{
				printf("%ld	%ld	%g	%d	%g\n", num, deno, ex, 0, gold);
				printf("%ld	%ld	%g	%d	%g\n", num, deno, ex, 1, gold);
				printf("\n");
			}
		}
		fflush(stdout);
	}
#endif
}
