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

#if 0
/* Evaluate the m,n'th eventually-golden polynomial at point x.
 */
double event_poly(unsigned long pfx, double x)
{
// printf("duuude n=%d x=%20.16g beta=\n", n, x, xn-acc);
	return 0;
}

/* Use midpoint bisection to find the single, unique
 * positive real zero of the n'th golden polynomial.
 */
double find_zero(unsigned long n, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = golden_poly(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}
#endif

int main(int argc, char* argv[])
{
	if (argc !=3)
	{
		fprintf(stderr, "Usage: %s <pfx> <cyc> <cyclen>\n", argv[0]);
		exit(1);
	}

	long pfx = atol(argv[1]);
	long cyc = atol(argv[2]);
	int cyclen = atoi(argv[3]);

	short cof[MAXCOF];
	get_coeffs(cof, pfx, cyc, cyclen);
}
