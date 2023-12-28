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
 *    (((moves & 1UL<<len) << 1) +1) / (1UL << (len-1))
 */
unsigned long good_index_map(unsigned long moves)
{
	unsigned long idx = 0;
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

void moves_to_rational(unsigned long moves, unsigned long* p, unsigned long* q)
{
	int nmov = bitlen(moves);
	unsigned long mask = 1UL << nmov;

	moves &= mask;
	moves <<= 1;
	moves |= 1UL;
	*p = moves;
	*q = mask>>1;
}

// Given a finite-orbit index, return the matching tree-walk (of left-
// right moves) to get to that index. Inverse of good_index_map()
unsigned long idx_to_moves(unsigned long idx)
{
	int len = 0;
	long bitseq = 0;
	// printf("enter idx=%ld\n", idx);
	while (NEG_ONE != idx)
	{
		bitseq |= 1UL << len;
		len++;
		while (0 == idx%2 && valid_gold_index(idx >>1))
		{
			idx >>= 1;
			len++;
		}

		long left = bracket_gold_left(idx);
		// printf("idx-to-moves: %ld |=> leader idx=%ld len=%d\n", left, idx, len);
		idx = left;
	}
	bitseq |= 1UL << len;
	return bitseq;
}

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

/* --------------------------- END OF LIFE ------------------------- */
