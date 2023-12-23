/*
 * bracket-finite.c
 *
 * Explore relationship between bracketing and mid-point moves.
 *
 * December 2023
 */

#include "selfie.c"

unsigned long move_gold_index(unsigned long moves, unsigned long idx)
{
	int nmov = bitlen(moves);

	for (int i=0; i< nmov; i++)
	{
		if (moves & 1UL)
			idx = move_gold_right(idx);
		else
			idx = move_gold_left(idx);

		moves >>= 1;
	}
	return idx;
}

int main(int argc, char* argv[])
{
	if (3 != argc) {
		fprintf(stderr, "Usage: %s <move> <idx>\n", argv[0]);
		exit(1);
	}
	long moves = atol(argv[1]);
	long idx = atol(argv[2]);

	print_gold_info(idx);
	long down = move_gold_index(moves, idx);
	print_gold_info(down);
}
