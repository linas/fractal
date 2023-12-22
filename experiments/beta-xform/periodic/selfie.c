/*
 * selfie.c
 *
 * Redesign of polynomial encoding and root finding.
 *
 * Notes about bit-encodings.
 * Orbits and dyadic fractions will use a little-endian encoding.
 * This is the same as 2-adic encoding, so that infinite strings
 * extend to the left.
 *
 * This means:
 *   -- bit-zero is the right-most bit.
 *   -- bit-one is <<1 to the left of bit zero.
 *   -- finite-length strings need an explicit length.
 *   -- printing will be left-to-right, so only the storage of
 *      the string is 2-adic, but not the printouts.
 *
 * This is a redesign and sometimes copy of code in irred-gold.c
 *
 * December 2023
 */

#define WORDLEN 64

/* Return length of bitstring. Same as ceil(log2(bitstr)). */
int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

// Print bitstring, so that bit zero appears left-most.
// Roughly speaking, this reverses the bit-order.
void print_dyadic(unsigned long bitseq, int len, char* pre, char* suf)
{
	printf("%s", pre);
	for (int i=0; i<len; i++)
	{
		printf("%ld", bitseq & 1UL);
		bitseq >>= 1;
	}
	printf (" \\\\ %d", len);
	printf("%s", suf);
}

/*
 * Return dyadic string corresponding to beta. This is obtained by
 * midpoint iteration, storing the bits in dyadic order, so that
 * first iteration is right-most bit. Note that right-most bit will
 * always be one, since the mid-point is always greater than 1/2.
 *
 * Finite orbits will return a finite string. Chaotic orbits might
 * use all 64 bits in the unsigned long.
 *
 * "Valid" integer betas will return the bitstring for (2n+1) but
 * in reversed order.
 */
unsigned long beta_to_dyadic(double beta)
{
	unsigned long bitseq = 0;
	double mid = 0.5*beta;
	for (int i=0; i < WORDLEN; i++)
	{
		if (0.5 <= mid)
		{
			bitseq |= 1UL << i;

#define MIDEPSI 1.0e-15
			// Apply rounding pressure, so as to favor finite iterates
			// over periodic ones. MIDEPSI=1e-15 seems to work at low
			// orders e.g. beta from index=6.
			mid -= 0.5-MIDEPSI;
		}
		mid *= beta;
	}
	return bitseq;
}

/* ================================================================= */

/* Implement the n'th golden polynomial. Return result from evaluating it.
 *
 * Polynomial is constructed from the bit string of (2n+1). Construction
 * is the bit-shift construction: the lowest powers of x are given by
 * right-most bits; highest powers are the left-most bits.
 */
double golden_poly(unsigned long idx, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	unsigned long bitstr = 2*idx+1;
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
static double find_zero(unsigned long idx, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = golden_poly(idx, mid);
	if (0.0 < fmid) return find_zero(idx, lo, mid);
	return find_zero(idx, mid, hi);
}

/* Return the beta value corresponding to the n'th golden polynomial.
 * Polynomial is constructed from the bit string of (2n+1).
 */
double golden_beta(unsigned long idx)
{
	return find_zero(idx, 1.0, 2.0);
}

/* ================================================================= */
/* Return true if `idx` is a self-describing index. This means that
 * the root of the golden_beta polynomial, when iterated by mid-point
 * iteration, reproduces the a dyadic bitstring, that, when bit-reversed,
 * gives back the index.
 *
 * The length of the orbit will be log_2(2n+1).
 *
 * This is returns exactly the same values as the theta() function,
 * except that it does not use a recursive algo to do it's work.
 */
bool valid_gold_index(unsigned long idx)
{
	double gold = golden_beta(idx);
	unsigned long dyad = beta_to_dyadic(gold);
	// printf("Index=%ld gold=%20.16g\n", idx, gold);
	// print_dyadic(dyad, 40, "Dyadic orbit=", "\n");

	unsigned long tno = 2*idx+1;
	int len = bitlen(tno);
	for (int i=0; i< len; i++)
	{
		if ((tno &1UL) != ((dyad>>(len-i-1)) &1UL)) return false;
		tno >>= 1;
	}
	return true;
}

/* ================================================================= */
// Return the dyadic that approximates the rational p/q,
// truncated to len bits. Rationals have infinite periodic bitseqs,
// so this just expands that bitseq, and then truncates it.
//
// Bit strings are stored right-to-left, so that the decimal point at
// far right. This is a 2-adic storage convention. For example, 1/2 and
// 1/4 and 1/8 are all stored as 0x1 but have length 1,2,3.
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

/* Return greatest common divisor of p and q */
unsigned long gcd(unsigned long p, unsigned long q)
{
	/* Euclid's algorithm for obtaining the gcd */
	while (0 != q)
	{
		unsigned long t = p % q;
		p = q;
		q = t;
	}
	return p;
}

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

	unsigned long cyc = 1UL;
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
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = event_poly(cof, mid);
	if (0.0 < fmid) return find_event_zero(cof, lo, mid);
	return find_event_zero(cof, mid, hi);
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
	unsigned long drat = rational_to_dyadic(p, q, WORDLEN);

	int badbit = compare_dyadics(dyad, drat);

	// For now, 40 bits seems OK-ish for the accuracy we've got on hand.
	if (0 <= badbit && badbit < 40) return false;
	return true;
}

/* --------------------------- END OF LIFE ------------------------- */
