/*
 * selfie-rational.c
 *
 * Code for the ultimately-periodic orbits. Stuff for handling rationals.
 * See selfie.c for groundwork basics.
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================= */
// Return the dyadic that approximates the rational p/q,
// truncated to len bits. Rationals have infinite periodic bitseqs,
// so this just expands that bitseq, and then truncates it.
//
// Bit strings are stored right-to-left, so that the decimal point at
// far right. This is a 2-adic storage convention. For example, 1/2 and
// 1/4 and 1/8 are all stored as 0x1 but have length 1,2,3.
//
// XXX this is junk and maybe should be deleted ??
//
unsigned long rational_to_dyadic(unsigned long p, unsigned long q, int len)
{
	unsigned long bitseq = 0;
	for (int i=0; i<len; i++)
	{
		p <<= 1;
		if (q <= p)
		{
			bitseq |= (1UL << i);
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
void get_event_cycle (unsigned long p, unsigned long q,
                      unsigned long *pfxp, unsigned long *cycp,
                      int *cyclenp)
{
	// printf("Enter get_event_cycle p/q = %ld/%ld\n", p,q);
	if (0 == p || 0 == q)
	{
		printf("Error: can't handle zero\n");
		*pfxp = 0; *cycp = 0; *cyclenp = 0;
		return;
	}

	unsigned long gcf = gcd(p+q, 2*q);
	unsigned long a = (p+q) / gcf;
	unsigned long b = (2*q) / gcf;
	// printf("Reduced a/b = %ld/%ld\n", a,b);

	int ell = 0;
	unsigned long br = b;
	while (br%2 == 0) { br >>= 1; ell++; }
	// printf("L = %d breduce = %ld\n", ell, br);

	int en = 1;
	while ((((1UL<<en) - 1UL) % br != 0) && en < WORDLEN) en++;

	if (WORDLEN == en)
	{
		printf("Error: cycle length overflow\n");
		*pfxp = 0; *cycp = 0; *cyclenp = 0;
		return;
	}
	// printf("N = cyclen = %d\n", en);

	unsigned long hi = (1UL<<en) - 1UL;
	unsigned long red = a * hi / br;
	// printf("reduce = %ld\n", red);

	unsigned long cyc = 0;
	while (((red-cyc) % hi != 0) && cyc < hi) cyc++;

	// printf("Cycle = cyc = %ld\n", cyc);
	unsigned long pfx = (red-cyc) / hi;

	if ((1UL<<ell) <= pfx)
	{
		printf("Error: Internal error of some kind\n");
		*pfxp = 0; *cycp = 0; *cyclenp = 0;
		return;
	}
	// printf("Prefix = pfx = %ld\n", pfx);
	*pfxp = pfx;
	*cycp = cyc;
	*cyclenp = en;
}

/* Max allowed coeficients for the eventually-periodic array. */
#define MAX_EVENT_COF 120

/* Convert (prefix,cyclic) coding pair to a list coefficients.
 * The bitseq for pfx must always start with one, so it's length is
 * unambiguous. The cyclic part is not bound by this constraint, so
 * it's length is explicitly specified.
 */
void get_event_coeffs(short *cof, long pfx, long cyc, int cyclen)
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

void print_event_coeffs(const short* cof, const char* pre, const char* suf)
{
	printf("%s", pre);
	for (int i=0; -100 < cof[i]; i++)
	{
		printf("%d ", cof[i]);
	}
	printf("%s", suf);
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
static double find_event_zero(short* cof, double lo, double hi)
{
	int nextra = 0;
	while (1)
	{
		double mid = 0.5 * (lo+hi);
		if (1.0e-15 > hi-lo) nextra++;
		if (3 < nextra) return mid; // run the loop three more times.
		double fmid = event_poly(cof, mid);
		if (0.0 == fmid) return mid;
		if (0.0 < fmid)
			hi = mid;
		else
			lo = mid;
	}
}

/* Utility wrapper for above. Get the polynomial zero for the
 * indicated rational p/q
 */
double event_gold(unsigned long p, unsigned long q)
{
	unsigned long pfx;
	unsigned long cyc;
	int cyclen;

	get_event_cycle (p, q, &pfx, &cyc, &cyclen);
	short cof[MAX_EVENT_COF];
	get_event_coeffs(cof, pfx, cyc, cyclen);
	double gold = find_event_zero(cof, 1.0, 2.0);
	return gold;
}

// ---------------------------------------------------------------------

/* Return index of the first miscompare. Return -1 if everything is OK. */
int compare_dyadics(unsigned long dya, unsigned long dyb)
{
	for (int i=0; i<WORDLEN; i++)
		if (((dya>>i) & 1UL) != ((dyb>>i) & 1UL))
			return i;

	return -1;
}

/*
 * Return true if the polynomial for rational p/q has a zero that,
 * when iterated, reproduces the bit-sequence for p/q (as a dyadic string).
 */
bool is_event_ok(unsigned long p, unsigned long q)
{
	double gold = event_gold(p, q);
	unsigned long dyad = beta_to_dyadic(gold);
	unsigned long drat = rational_to_dyadic(p+q, 2*q, WORDLEN);

	int badbit = compare_dyadics(dyad, drat);

	// For now, 40 bits seems OK-ish for the accuracy we've got on hand.
	if (0 <= badbit && badbit < 40) return false;
	return true;
}

/* --------------------------- END OF LIFE ------------------------- */
