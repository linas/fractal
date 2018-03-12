/*
 * tribo.c
 *
 * Generalized Tribonacci numbers -- these are the same as
 * the polynomial roots.
 *
 * Linas Vepstas February 2018
 */
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	int nmax = 50;
	unsigned long seq[nmax];

	if (argc<3)
	{
		fprintf(stderr, "Usage: %s len bits...\n", argv[0]);
		fprintf(stderr, "Example Narayana: %s 3 1 0 1\n", argv[0]);
		fprintf(stderr, "Example Tribonacci: %s 3 1 1 1\n", argv[0]);
		exit(1);
	}

	int ilen = atoi(argv[1]);
	unsigned long mask[ilen];
	for (int j=0; j<ilen; j++)
	{
		seq[j] = 0;
		mask[j] = atoi(argv[2+j]);
		printf("# Mask: %d	%lu\n", j, mask[j]);
	}
	seq[ilen-1] = 1;

	for (int i=ilen; i< nmax; i++)
	{
		unsigned long sum = 0;
		for (int j=0; j<ilen; j++) sum += mask[j] * seq[i-j-1];

		seq[i] = sum;

		double rat = ((double) seq[i]) / ((double) seq[i-1]);
		printf("%d	%lu	%20.18g\n", i, seq[i], rat);
	}
}
