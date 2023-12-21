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
	printf (" \\%d", len);
	printf("%s", suf);
}

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
		bitseq <<= 1;
		p <<= 1;
		if (q <= p)
		{
			bitseq |= 1UL;
			p -= q;
		}
	}
	return bitseq;
}

/* --------------------------- END OF LIFE ------------------------- */
