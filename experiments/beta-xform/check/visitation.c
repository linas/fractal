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
double midp[DEPTH];

void iter(double beta)
{
	double midpnt = 0.5*beta;
	double obn = 1.0;
	for (int i=0; i<DEPTH; i++)
	{
		midp[i] = midpnt;

		// Search
		double least = 42;
		int rig = -1;
		for (int j=0; j<i; j++)
		{
			if (least < midp[j]) continue;
			if (midpnt < midp[j] && midp[j] < least)
			{
				least = midp[j];
				rig = j;
			}
		}
printf("i=%d mp=%g least=%g rig=%d\n", i, midpnt, least, rig);

		// Step to next midpoint
		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		if (obn < 1e-15) break;
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

	// Setup
	for (int i=0; i<DEPTH; i++) visit[i] = 0;

	iter(beta);
}
