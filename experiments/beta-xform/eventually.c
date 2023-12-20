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
 * Polynomial coefficients already decoded in cof.
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

/* Make sure that midpoint iteration gives same bitstring
 * as the encoding w/ pfx, cyc
 */
bool validate_orbit(double beta, long pfx, long cyc, int cyclen)
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
			printf("Error out at pfx %d %d %d\n", i, pfbit, bit);
			return false;
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
			printf("Error out at cyclic %d %d %d\n", i, cybit, bit);
			return false;
		}
	}
	return true;
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <pfx> <cyc> <cyclen>\n", argv[0]);
		exit(1);
	}

	long pfx = atol(argv[1]);
	long cyc = atol(argv[2]);
	int cyclen = atoi(argv[3]);

	short cof[MAXCOF];
	get_coeffs(cof, pfx, cyc, cyclen);
	print_coeffs(cof);
	double gold = find_ezero(cof, 1.0, 2.0);

	printf("found %g\n", gold);
	validate_orbit(gold, pfx, cyc, cyclen);
}
