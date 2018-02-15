/*
 * tribo.c
 *
 * Tribonacci foolishness
 * Linas Vepstas February 2018
 */
#include <stdio.h>

int main ()
{
	int nmax = 50;
	int seq[nmax];

	seq[0] = 1;
	seq[1] = 0;
	seq[2] = 1;

	for (int i=3; i< nmax; i++)
	{
		int sum = 0;
		for (int j=0; j<3; j++) sum += seq[i-j-1];
		seq[i] = sum;

		double rat = ((double) seq[i]) / ((double) seq[i-1]);
		printf("%d	%d	%g\n", i, seq[i], rat);
	}
}
