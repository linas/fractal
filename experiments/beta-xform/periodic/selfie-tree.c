/*
 * selfie-tree.c
 *
 * Code for tree moves. This includes the code for tie "good index map".
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================= */
/*
 * The "good index map". Given a sequence of L,R moves on the full
 * binary tree, return the corresponding valid finite index. This is
 * the "G" function in the paper.
 *
 * Return the front-center index from the bracketing sequence. The
 * bracketing sequence is a map of valid indexes to the binary tree.
 * The encoding is with move_gold_left() and move_gold_right(), which
 * respectively either double the index (moving left), or find the next
 * leader (when moving right).
 *
 * The fronts can be understood to be encoded in the sequence number
 * itself, taking that number as series of left-right moves down the
 * tree. The moves are done by reading the bits in the binary rep of
 * the number, from left to right. The left-most bit is necessarily
 * a one, so the first move is necessarily R, and so we start with
 * index=0, which is beta=1, which is the left feeder to the root.
 * after that, we are at the root, and so everything afterwards is just
 * an encoding of all-possible left-right moves.
 *
 * Thus:
 * root == 1/2 == binary 1 length=1
 * left move == L == 1/4 == binary 10 length = 2
 * right move == R == 3/4 == binary 11 length = 2
 * LL == 1/8 == 100, length=3
 * LR == 3/8 == 101, length=3
 * RL == 5/8 == 110, length=3
 * RR == 7/8 == 111, length=3
 *
 * So this is a "one-leading string", which is the opposite of the
 * dyadic encoding used in other parts of this file. The nice thing
 * about the one-leading encoding is that an explicit length is not
 * needed; just lop off the MSB bit, and the rest of the string follows.
 *
 * The matching dyadic fraction is obtained by trashing the leading
 * bit, then (2n+1), then divide by 2^len, as follows:
 *    (((moves & 1UL<<len) << 1) +1) / (1UL << len)
 *
 * inverse given by idx_to_moves()
 */
unsigned long good_index_map(unsigned long moves)
{
	// Special-case the endpoints
	// No moves is far left beta=1
	// All-right moves is RRRR=-1 which overflows but should be beta=2
	if (0UL == moves)
		return NEG_ONE;
	if (NEG_ONE == moves)
		return 0UL;

	unsigned long idx = 0UL;
	int nmov = bitlen(moves);
	for (int i=0; i< nmov; i++)
	{
		if (moves & (1UL<<(nmov-1-i)))
			idx = move_gold_right(idx);
		else
			idx = move_gold_left(idx);

		if (MAXIDX < idx) return -2;
	}
	return idx;
}

// Given a finite-orbit index, return the matching tree-walk (of left-
// right moves) to get to that index. Inverse of good_index_map()
unsigned long idx_to_moves(unsigned long idx)
{
	// Far left is beta=1 which is no moves at all, so zero.
	// Far right is beta=2 which index 0, which is all right moves, which is -1
	if (NEG_ONE == idx) return 0;
	if (0 == idx) return -1;
	if (MAXIDX < idx) return -2;

	int len = 0;
	long bitseq = 0;
	// printf("idx-to-moves: --- Enter idx=%ld\n", idx);
	while (NEG_ONE != idx)
	{
		while (0 == idx%2 && valid_gold_index(idx >>1))
		{
			idx >>= 1;
			len++;
		}
		bitseq |= 1UL << len;
		len++;

		long left = bracket_gold_left(idx);
		// printf("idx-to-moves: %ld |=> leader idx=%ld len=%d\n", left, idx, len);
		idx = left;
	}
	return bitseq;
}

/* ================================================================= */

// Print moves on the full binary tree.
void print_moves(unsigned long moves, char* pre, char* suf)
{
	int len = bitlen(moves);
	printf("%s", pre);
	for (int i=0; i<len; i++)
	{
		int bit = (moves >> (len-i-1)) & 0x1L;
		if (bit) printf("R");
		else printf("L");
	}
	printf (" /%d", len);
	printf("%s", suf);
}

/* ================================================================= */

// Given a tree location encoded with canonical integer, return the corresponding
// rational for that location. The root of the tree is 1/2, the fisrt row down
// is 1/4 and 3/4, then next row is 1/8, 3/8, 5/8, 7/8, and so on.
// This is the canonical dyadic binary tree.
// The inverse function is rational_to_moves(), below.
void moves_to_rational(unsigned long moves, unsigned long* p, unsigned long* q)
{
	// No moves is far left, beta=1
	if (0 == moves)
	{
		*p = 0;
		*q = 1;
		return;
	}

	int nmov = bitlen(moves);
	unsigned long deno = 1UL << nmov;
	unsigned long mask = deno - 1UL;
	mask >>= 1;

	moves = moves & mask;
	moves <<= 1;
	moves |= 1UL;
	*p = moves;
	*q = deno;
}

// Return true, if q is a pure power of two.
bool is_power_of_two(unsigned long q)
{
	if (1 == q) return true;
	while (1 < q && 0 == q%2) q >>= 1;
	return 1 == q;
}

// Return the bit-sequence that approximates the rational p/q,
// truncated to len bits. Rationals have infinite periodic bitseqs,
// so this just expands that bitseq, and then truncates it.
//
// Bit strings are read left-to-right, with leading one taking
// indicating location of decimal point. This gives the canonical
// mapping for dyadic rationals as integers.
// For example, 1/2 len=1 -> 1
// For example, 1/4 len=2 -> 10 and 3/4 = 11
// and  3/8 len=3 -> 101
//
// A note about dyadic ratinoals and length. If q is a power of two,
// then the iteration will reach zero, and L moves will be appended
// until `len` is reached. This may or may not be what you want. If
// it is what you want, then fine. Otherwise, "don't do that". Specify
// a length that is appropriate. The following will do:
//
//    if (is_power_of_two(q) len = bitlen(q) - 2;
//
// This is the inverse of move_to_rational(), above.
unsigned long rational_to_moves(unsigned long p, unsigned long q, int len)
{
	unsigned long gcf = gcd(p, q);
	p /= gcf;
	q /= gcf;

	if (0 == len) return 1UL;
	if (-1 == len) return 0UL;

	unsigned long bitseq = 1UL;
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
