/*
 * visitation.c
 *
 * Visitation function, from the diary.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DEPTH 5000
unsigned long visit[DEPTH];
unsigned long rigb[DEPTH];
double midp[DEPTH];
int sorted[DEPTH];

int cmp(const void* ida, const void* idb)
{
	int ia = *((int*) ida);
	int ib = *((int*) idb);
	return (midp[ia] < midp[ib]) ? -1: 1;
}

#define WORDLEN 64
#define MAXIDX (1UL<<(WORDLEN-1))
#define NEG_ONE ((unsigned long)(-1L))

void iter(double beta)
{
	// Setup
	for (int i=0; i<DEPTH; i++) visit[i] = NEG_ONE;
	visit[0] = 0;
	rigb[0] = 1;

	double midpnt = 0.5*beta;
	midp[0] = midpnt;
	double obn = 1.0;
	int tot = 0;
	for (int i=1; i<DEPTH; i++)
	{
		// Step to next midpoint
		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
#define EPS 1e-15
		// if (obn < EPS) break;

		midp[i] = midpnt;

		// Search
		double rmost = 4200;
		int rig = -1;
		for (int j=0; j<i; j++)
		{
			if (rmost < midp[j]) continue;
			if (midpnt < midp[j] && midp[j] < rmost)
			{
				rmost = midp[j];
				rig = j;
			}
		}

		if (MAXIDX < rigb[rig]) break;

		// Record right-hand boundary, and update it.
		visit[i] = rigb[rig];
		rigb[i] = 2*rigb[rig];
		rigb[rig] = 2*rigb[rig]+1;
		// printf("i=%d mp=%g least=%g rig=%d -- visit= %ld\n",
		//	i, midpnt, rmost, rig, visit[i]);
		tot = i;
	}
	tot++;

	// Sort into sequential order.
	for (int i=1; i<DEPTH; i++) sorted[i] = i;

	qsort(sorted, tot, sizeof(int), cmp);
}

/* Return length of bitstring. Same as ceil(log2(bitstr)). */
int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/* Return canonical dyadic associated with canonical integer
 * in the canonical tree numbering.The map is
 *   1 |--> 1/2
 *   2 |--> 1/4
 *   3 |--> 3/4
 *   4 |--> 1/8
 *   5 |--> 3/8
 * etc.
 */
double canonical_dyadic(unsigned long n)
{
	int len = bitlen(n);
	unsigned long base = 1UL<<len;
	unsigned long nb = n - base/2;
	unsigned long nd = 2*nb+1;
	double dya = ((double) nd) / ((double) base);
	return dya;
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);

	iter(beta);

	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 0.0, 0.0);
	for (int i=0; i<DEPTH; i++)
	{
		if (NEG_ONE == visit[i]) break;
		int j = sorted[i];
		unsigned long v = visit[j];
		double dya = canonical_dyadic(v);
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, midp[j]);
	}
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 1.0, 1.0);

#ifdef STEPPY
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 0.0, 0.0);
	double yprev = 0.0;
	for (int i=0; i<DEPTH; i++)
	{
		if (NEG_ONE == visit[i]) break;
		int j = sorted[i];
		unsigned long v = visit[j];
		double dya = canonical_dyadic(v);
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, yprev);
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, midp[j]);
		yprev = midp[j];
	}
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 1.0, 1.0);
#endif
}
