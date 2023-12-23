/*
 * bracket-finite.c
 *
 * Explore relationship between bracketing and mid-point moves.
 *
 * December 2023
 */

#include "selfie.c"

/*
 * Implement the result of making a sequence of left-right moves down
 * the binary tree, starting from the root.
 */
unsigned long move_gold_index(unsigned long moves)
{
	unsigned long idx = 0;
	int nmov = bitlen(moves);
	for (int i=0; i< nmov; i++)
	{
		if (moves & (1UL<<(nmov-1-i)))
			idx = move_gold_right(idx);
		else
			idx = move_gold_left(idx);
	}
	return idx;
}

int main(int argc, char* argv[])
{
#ifdef MANUAL_EXPLORER
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <move>\n", argv[0]);
		exit(1);
	}
	long moves = atol(argv[1]);

	printf("Moves from index zero:\n");
	long down = move_gold_index(moves);
	print_gold_info(down);

	printf("Just the index:\n");
	print_gold_info(moves);
#endif

	// Unit test
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	long maxind = 1UL << maxord;

	for (long mov=0; mov< maxind; mov++)
	{
		bool ok = valid_gold_index(mov);
		long idx = move_gold_index(mov);

		printf("%ld	%d	%ld\n", mov, ok, idx);

		// long lead = gold_leader(mov);
	}
}
