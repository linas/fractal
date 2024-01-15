/*
 * deriv.c
 *
 * Explorederivative of polynomials.
 *
 * January 2024
 */

#include "selfie.c"

int main(int argc, char* argv[])
{
#define MANUAL_EXPLORER
#ifdef MANUAL_EXPLORER
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxorder>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	long maxidx = 1UL << maxord;

	for (long idx = 1; idx < maxidx; idx++)
	{
		bool ok = valid_gold_index(idx);
		if (!ok) continue;

		double gold = golden_beta(idx);
		printf("Its %ld %g\n", idx, gold);
	}
#endif
}
