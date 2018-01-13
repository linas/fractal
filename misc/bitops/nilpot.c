
/*
 * nilpot.c
 *
 * Explore the niplotent irrationals.
 * January 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s NMAX\n", argv[0]);
		exit (1);
	}
	int nmax = atoi(argv[1]);

#define NBINS 1501
	int cnt[NBINS];
	for (int i=0; i< NBINS; i++) cnt[i] = 0;

	for (int n=2; n< nmax; n++)
	{
	}
}
