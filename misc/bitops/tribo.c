/*
 * tribo.c
 *
 * Tribonacci foolishness
 * Linas Vepstas February 2018
 */
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	int nmax = 50;
	unsigned long seq[nmax];
	seq[0] = 1;
	seq[1] = 0;
	seq[2] = 1;

	if (argc<3)
	{
		fprintf(stderr, "Usage: %s len vals...", argv[0]);
		exit(1);
	}

	int ilen = atoi(argv[1]);
	for (int j=0; j<ilen; j++)
	{
		seq[j] = atoi(argv[2+j]);
		printf("# Initial sequence: %d	%lu\n", j, seq[j]);
	}

	for (int i=ilen; i< nmax; i++)
	{
		unsigned long sum = 0;
		for (int j=0; j<ilen; j++) sum += seq[i-j-1];
		seq[i] = sum;

		double rat = ((double) seq[i]) / ((double) seq[i-1]);
		printf("%d	%lu	%20.18g\n", i, seq[i], rat);
	}
}
