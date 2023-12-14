/*
 * irred-allowed.c
 *
 * Verify the properties of the allowed sequence.
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "irred-gold.c"

// =================================================================

int main(int argc, char* argv[])
{
	// Obtain index from command line.
	if (1 == argc) {
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}

	long maxidx = atol(argv[1]);

	malloc_gold(maxidx);
	malloc_index_cache(maxidx);

	// This should always be silent.
	for (long i = 0; i< maxidx/2; i++)
	{
		bool valid = is_valid_index(i);
		if (false == valid) continue;

		bool v2 = is_valid_index(2*i);
		if (false == v2)
		{
			printf("Error: bad doubling %ld\n", i);
		}

		if (1 == i%2)
		{
			bool v3 = is_valid_index((i-1)/2);
			if (false == v3)
			{
				printf("Error: bad tripling %ld\n", i);
			}
		}
	}
}
