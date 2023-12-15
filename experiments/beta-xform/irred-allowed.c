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

	long maxord = atol(argv[1]);
	long maxidx = 1UL << maxord;

	malloc_gold(4*maxidx);
	malloc_index_cache(4*maxidx);

	// This should always be silent.
	for (long i = 0; i< maxidx; i++)
	{
		bool valid = is_valid_index(i);
		bool stop = is_stopper(i);
		if (stop == valid)
			printf("Error: bad stop at %ld valid=%d stop=%d\n", i, valid, stop);

		if (false == valid)
		{
			bool vodd = is_valid_index(2*i+1);
			if (vodd)
				printf("Bad odd i=%ld\n", i);
			for (long j=1; ; j++)
			{
				long k = 1UL << j;
				k *= 2*i+1;
				if (maxidx <= k) break;
				bool prod = is_valid_index(k);
				if (prod)
					printf("Huh prod odd i=%ld j=%ld\n", i,j);
			}
		}

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

	// print_stoppers(256);
	print_stoppers(2048);

#if 0
	for (int p=0; p<20; p++)
	{
		long b = 1L << p;
		for (long w = 1; w< p; w++)
		{
			long d = b - (2*w-1);
			bool vd = is_valid_index(d);
			if (false == vd)
			{
				printf("fail at p=%d base=%ld w=%ld del=%ld d=%ld\n", p, b, w, 2*w-1, d);
			}
		}
	}
#endif

#if 0
	// Sum of Moreau
	printf("yo %d is %ld\n", 1, valid_index_cache(1));
	printf("yo %d is %ld\n", 3, valid_index_cache(3));
	printf("yo %d is %ld\n", 6, valid_index_cache(6));
	printf("yo %d is %ld\n", 12, valid_index_cache(12));
	printf("yo %d is %ld\n", 21, valid_index_cache(21));
	printf("yo %d is %ld\n", 39, valid_index_cache(39));
	printf("yo %d is %ld\n", 69, valid_index_cache(69));
	printf("yo %d is %ld\n", 125, valid_index_cache(125));
#endif
}
