/*
 * good-rational.c
 *
 * Create a graph of the bracket tree.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-rational.c"
#include "selfie-tree.c"

// cheap hack
void low_guess(double rat, int* p, int* q)
{
	for (int i=2; i<16386; i++)
	{
		if ((fabs(fmod(rat*i, 1.0)) < 1e-6) ||
		    (fabs(1.0-fmod(rat*i, 1.0)) < 1e-6))
		{
			*q = i;
			*p = floor(rat*i + 1e-5);
			return;
		}
	}
}

int main(int argc, char* argv[])
{
#define RATIONAL_EXPLORE
#ifdef RATIONAL_EXPLORE
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <p> <q> <max-len>\n", argv[0]);
		exit(1);
	}
	int p = atoi(argv[1]);
	int q = atoi(argv[2]);
	int len = atoi(argv[3]);
	printf("\n");

	double rat = ((double) p) / ((double) q);
	printf("Input rational = %d/%d len=%d as float=%18.16g\n", p, q, len, rat);

	long moves = rational_to_moves(p, q, len);
	print_moves(moves, "Tree moves=(", ")\n");

	unsigned long pp, qq;
	moves_to_rational(moves, &pp, &qq);
	double rrat = ((double) pp) / ((double) qq);
	printf("Reconstructed rational = %ld/%ld = %18.16g\n", pp, qq, rrat);

	unsigned long idx = good_index_map(moves);
	double gold = golden_beta(idx);
	printf("Bracket index %ld  beta=%20.16g\n", idx, gold);

	unsigned long oradic = beta_to_dyadic(gold);
	print_dyadic(oradic, 63, "Actual orbit: ", "\n");

	unsigned long oridx = beta_to_index(gold);
	printf("Reconstructed index from orbit: %ld\n", oridx);
	if (oridx != idx)
	{
		double rgold = golden_beta(oridx);
		printf("Warning: indexes don't match! diff=%g\n", gold-rgold);
	}

	unsigned long tno = 2UL * idx + 1UL;
	moves_to_rational(tno, &pp, &qq);
	double orat = ((double) pp) / ((double) qq);
	printf("Orbit rational = %ld/%ld = %18.16g\n", pp, qq, orat);

	int po=0, qo=0;
	low_guess(orat, &po, &qo);
	printf("guess its orbit is %d/%d\n", po,qo);
	printf("\n");

#endif
}

/* --------------------------- END OF LIFE ------------------------- */
