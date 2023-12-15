/*
 * irred-height.c
 *
 * Exploration of the height of the leaders.
 *
 * December 2023
 */

#include "irred-gold.c"

int main(int argc, char* argv[])
{
	int maxord = atoi(argv[1]);

	long nmax = 1 << maxord;
	malloc_gold(nmax);

	for (long idx=0; idx<nmax; idx++)
	{
		if (false == is_valid_index(idx)) continue;
		int height = 0;
		long next = 2*idx+1;
		do
		{
			if (is_valid_index(next))
			{
				double gold = find_gold(idx);
				double lead = find_gold(next);
				printf("%g	%g	%ld	%d\n", gold, lead, idx, height);
				break;
			}
			height ++;
			next <<= 1;
		}
		while (next <nmax);
	}
}

/* --------------------------- END OF LIFE ------------------------- */
