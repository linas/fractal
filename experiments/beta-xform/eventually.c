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

#include "irred-gold.c"

/* Max allowed coeficients */
#define MAXCOF 120

int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/* Convert (prefix,cyclic) coding pair to a list coefficients.
 * The bitseq for pfx must always start with one, so it's length is
 * unambiguous. The cyclic part is not bound by this constraint, so
 * it's length is explicitly specified.
 */
void get_coeffs(short *cof, long pfx, long cyc, int cyclen)
{
	int pfxlen = bitlen(pfx);
	for (int i=0; i< pfxlen; i++)
		cof[i] = pfx >> (pfxlen-i-1) & 1UL;

	for (int i=0; i< cyclen; i++)
		cof[i+pfxlen] = cyc >> (cyclen-i-1) & 1UL;

	cof[pfxlen+cyclen] = -666; // end of string marker;

	cof[cyclen-1] += 1;
	for (int i=0; i< pfxlen; i++)
		cof[cyclen+i] -= pfx >> (pfxlen-i-1) & 1UL;
}

void print_coeffs(short* cof)
{
	for (int i=0; -100 < cof[i]; i++)
	{
		printf("%d ", cof[i]);
	}
	printf("\n");
}

/*
 * Evaluate the eventually-golden polynomial at point x.
 * Polynomial coefficients must be provided in cof.
 */
double event_poly(short* cof, double x)
{
	int clen = 0;
	for (int i=0; -100 < cof[i]; i++) clen++;

	double f = 0.0;
	double xn = 1.0;
	for (int i=0; i<clen; i++)
	{
		f += cof[clen-i-1] * xn;
		xn *= x;
	}
	f = xn - f;
// printf("duuude x=%20.16g beta=\n", x, f);
	return f;
}

/* Use midpoint bisection to find the single, unique
 * positive real zero of the polynomial.
 */
double find_ezero(short* cof, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = event_poly(cof, mid);
	if (0.0 < fmid) return find_ezero(cof, lo, mid);
	return find_ezero(cof, mid, hi);
}

/* Utility wrapper for above. Get the polynomial zero for the
 * indicated prefix and cycle
 */
double event_gold(long pfx, long cyc, int cyclen)
{
	short cof[MAXCOF];
	get_coeffs(cof, pfx, cyc, cyclen);
	double gold = find_ezero(cof, 1.0, 2.0);
	return gold;
}

// ---------------------------------------------------------------------
/* Make sure that midpoint iteration gives same bitstring
 * as the encoding w/ pfx, cyc.
 * Return index of first disagreement.
 */
int validate_orbit(double beta, long pfx, long cyc, int cyclen, bool prt)
{
#define MAXBITS 60

	int pfxlen = bitlen(pfx);

	double mid = 0.5*beta;
	for (int i=0; i < pfxlen; i++)
	{
		int bit = 0;
		if (0.5 < mid) bit = 1;
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;

		int pfbit = pfx >> (pfxlen-i-1) & 1UL;
		if (bit != pfbit)
		{
			if (prt) printf("Error out at pfx %d %d %d\n", i, pfbit, bit);
			return i;
		}
	}
	for (int i=0; i < MAXBITS-pfxlen; i++)
	{
		int bit = 0;
		if (0.5 < mid) bit = 1;
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;

		int j = i%cyclen;
		int cybit = cyc >> (cyclen-j-1) & 1UL;
		if (bit != cybit)
		{
			if (prt) printf("Error out at cyclic %d %d %d\n", i, cybit, bit);
			return i+pfxlen;
		}
	}
	return 0;
}

/* Return true if the beta generates the expected orbit */
bool is_orbit_ok(long pfx, long cyc, int cyclen)
{
	double gold = event_gold(pfx, cyc, cyclen);

	int badbit = validate_orbit(gold, pfx, cyc, cyclen, false);
	if (0 < badbit && badbit < 40) return false;
	return true;
}

/* Iterate the midpoint and return the sequence */
long orbit_to_bitseq(double beta)
{
#define MAXBITS 60
	long bitseq = 0;
	double mid = 0.5*beta;
	for (int i=0; i < MAXBITS; i++)
	{
		bitseq <<= 1;
		if (0.5 < mid) bitseq |= 1UL;
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;
	}
	return bitseq;
}

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
	short cof[MAXCOF];
	get_coeffs(cof, pfx, cyc, cyclen);
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
	return is_valid_index(idx);
}

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

void print_bitseq(long bitseq, char* pre, char* suf)
{
	printf("%s", pre);
	for (int i=0; i< MAXBITS; i++)
		printf("%ld", (bitseq >> (MAXBITS-1-i)) & 1UL);

	printf("%s", suf);
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

// ---------------------------------------------------------------------

void print_debug_info(long pfx, long cyc, int cyclen)
{
	bool pok = is_prefix_ok(pfx, cyc, cyclen);
	if (false == pok)
		printf("Error: prefix is reducible\n");

	short cof[MAXCOF];
	get_coeffs(cof, pfx, cyc, cyclen);
	printf("Coeffs are: ");
	print_coeffs(cof);
	double gold = find_ezero(cof, 1.0, 2.0);

	printf("Prefix = %ld cycle=%ld len=%d\n", pfx, cyc, cyclen);
	double rat = event_to_double(pfx, cyc, cyclen);
	printf("Orbit code=%20.16g\n", rat);

	printf("Root %20.16g\n", gold);
	int badbit = validate_orbit(gold, pfx, cyc, cyclen, true);
	if (0 < badbit && badbit < 40)
		printf("Error: bad orbit at %d\n", badbit);

	long orbits = orbit_to_bitseq(gold);
	long events = event_to_bitseq(pfx, cyc, cyclen);
	print_bitseq(orbits, "orbit= ", "\n");
	print_bitseq(events, "event= ", "\n");

#define ONE_SIXTH
#ifdef ONE_SIXTH
	long convbits = orbit_to_bitseq(1.438430227468187);
	print_bitseq(convbits, "convb= ", "\n");
	long convbit2 = orbit_to_bitseq(1.438385389665382);
	print_bitseq(convbit2, "conv2= ", "\n");
#endif

#ifdef ONE_THIRD
	long convbits = orbit_to_bitseq(1.558993131521904);
	print_bitseq(convbits, "convb= ", "\n");
	long convbit2 = orbit_to_bitseq(1.558956158884443);
	print_bitseq(convbit2, "conv2= ", "\n");
#endif
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
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <pfx> <cyc> <cyclen>\n", argv[0]);
		exit(1);
	}

	long pfx = atol(argv[1]);
	long cyc = atol(argv[2]);
	int cyclen = atoi(argv[3]);

	print_debug_info(pfx, cyc, cyclen);
#endif

// #define BULK_LISTING
#ifdef BULK_LISTING
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s <maxpfx> <maxcyclen>\n", argv[0]);
		exit(1);
	}

	long maxpfx = atol(argv[1]);
	int maxcyclen = atoi(argv[2]);

#ifdef FAILED_HYPOTHESIS
	int maxord = bitlen(maxpfx) + maxcyclen;
	long maxidx = 1UL << maxord;
	malloc_gold(maxidx);
#endif

#define NROOTS 1000
	double roots[NROOTS];
	int totfnd = 0;
	for (int cyclen =2; cyclen <=maxcyclen; cyclen++)
	{
		long maxcyc = (1<<cyclen) - 1UL;
		for (long cyc = 1; cyc <maxcyc; cyc++)
		{
			for (long pfx = 1; pfx <maxpfx; pfx++)
			{
				if (false == is_prefix_ok(pfx, cyc, cyclen)) continue;

#ifdef FAILED_HYPOTHESIS
				// This did not work. Hope was that we could deduce stuff from
				// the finite strings, but no such luck.
				bool orbok = is_orbit_ok(pfx, cyc, cyclen);
				if (false == is_valid_finite(pfx, cyc, cyclen))
				{
					if (orbok) printf ("Yoooooooooo bad finite but OK orbit\n");
				}
				else
				{
					if (!orbok) printf ("Gooodddd finite but bad orbit!\n");
				}
#endif
				if (false == is_orbit_ok(pfx, cyc, cyclen)) continue;

				double gold = event_gold(pfx, cyc, cyclen);
				printf("Found (%ld, %ld/%d) = %20.16g\n", pfx, cyc, cyclen, gold);
				roots[totfnd] = gold;
				totfnd ++;
			}
			printf("------\n");
		}
		printf("=======\n");
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
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s <maxpfx> <maxcyclen>\n", argv[0]);
		exit(1);
	}

	long maxpfx = atol(argv[1]);
	int maxcyclen = atoi(argv[2]);

	for (int cyclen =2; cyclen <=maxcyclen; cyclen++)
	{
		long maxcyc = (1<<cyclen) - 1UL;
		for (long cyc = 1; cyc <maxcyc; cyc++)
		{
			for (long pfx = 1; pfx <maxpfx; pfx++)
			{
				if (false == is_prefix_ok(pfx, cyc, cyclen)) continue;
				bool orbok = is_orbit_ok(pfx, cyc, cyclen);

				double rat = orbit_to_double(pfx, cyc, cyclen);
				double gold = event_gold(pfx, cyc, cyclen);
				double good = 0.0;
				if (orbok) good = gold;
				double bad = 0.0;
				if (false==orbok) bad = gold;
				double fin = finite_to_double(pfx, cyc, cyclen);
				printf("%ld	%ld	%d	%g	%g	%g	%g\n", pfx, cyc, cyclen, rat, good, bad, fin);
			}
		}
	}
#endif
}
