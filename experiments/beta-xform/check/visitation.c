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

#define DEPTH 300
unsigned long visit[DEPTH];
unsigned long rigb[DEPTH];
double midp[DEPTH];

void iter(double beta)
{
	// Setup
	for (int i=0; i<DEPTH; i++) visit[i] = 0;
	visit[0] = 1;
	rigb[0] = 1;

	double midpnt = 0.5*beta;
	midp[0] = midpnt;
	double obn = 1.0;
	for (int i=1; i<DEPTH; i++)
	{
		// Step to next midpoint
		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		if (obn < 1e-15) break;

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

		visit[i] = rigb[rig];
		rigb[i] = 2*rigb[rig];
		rigb[rig] = 2*rigb[rig]+1;
		// printf("i=%d mp=%g least=%g rig=%d -- visit= %ld\n",
		//	i, midpnt, rmost, rig, visit[i]);
	}
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
}
