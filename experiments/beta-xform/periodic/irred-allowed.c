/*
 * irred-allowed.c
 *
 * Verify the properties of the allowed sequence.
 * Basically, this is a unit test for the code, making sure nothing
 * crazy is going on, and that the valid index sequence works as intended.
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
		fprintf(stderr, "Usage: %s <order>\n", argv[0]);
		exit(1);
	}

	long maxord = atol(argv[1]);
	long maxidx = 1UL << maxord;

	malloc_gold(4*maxidx);
	malloc_index_cache(4*maxidx);

#if 0
	long k = 1;
	long mp = 10;
	do
	{
		bool valid = is_valid_index(mp);
		bool vh = is_valid_index(mp/2);
		if (!valid || vh)
			printf("oh nooooo %ld\n", k);
		k++;
		long p = 1UL << k;
		mp = p * (2*p+1);
	} while (mp < 4*maxidx);
	printf("done up to k=%ld\n", k);
	exit(0);
#endif

#if 0
	// Handy dandy explorer for leaders of 2.
	long ldr = 2;
	long xpl = 2;
	while (0 < ldr)
	{
		printf("ldr = %ld  %ld\n", ldr, xpl);
		ldr = find_leader(ldr);
		xpl = 2UL * (2UL * xpl + 1UL);
	}
	exit(0);
#endif

#if 1
	// Print leaders

	int upcnt[32];
	for (int i=0; i<32; i++) upcnt[i] = 0;

	for (int nu = 2; nu< maxord; nu++)
	{
		int ntall = 0;
		long istart = (1UL<<(nu-2));
		long iend = 2*istart;
		for (long i = istart; i< iend; i++)
		{
			bool valid = is_valid_index(i);
			if (false == valid) continue;
			// long ldr = find_leader(i);
			long ldr = gold_leader(i);
			if (ldr <= 0) { printf ("Error: overflow at %ld\n", i); }

			long pheight = ldr / (2*i+1);
			int height = bitlen(pheight) - 1;
			if (0 < height) ntall ++;

			upcnt[nu+height] ++;

			// printf("ord = %d idx = %ld ldr = %ld height=%d\n", nu, i, ldr, height);
		}
		printf("ord = %d ntal = %d\n", nu, ntall);
	}

	printf("\nPopulation at each height:\n");
	for (int nu = 2; nu< maxord; nu++)
		printf("ord = %d pop = %d\n", nu, upcnt[nu]);

	printf("\n");
	for (int nu = 2; nu< maxord; nu++)
		printf("%d, ", upcnt[nu]);
	printf("\n");
	exit(0);
#endif

// #define SELF_CONSISTENCY_CHECK
#ifdef SELF_CONSISTENCY_CHECK
	// Except for the printf, this should always be silent.
	printf("Check self-consistency of the irred-gold.c codebase\n");
	for (long i = 0; i< maxidx; i++)
	{
		bool valid = is_valid_index(i);

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
	// print_stoppers(2048);
#endif

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
