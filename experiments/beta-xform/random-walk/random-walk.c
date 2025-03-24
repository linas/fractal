/*
 * random-walk.c
 *
 * Gaussian random walk.
 * Linas Vepstas 24 March 2025
 */

#include <stdlib.h>
#include <stdio.h>

/* Number of histogram bins. */
#define HISTO_NBINS 302

/* Number of samples */
#defin NSAMP 5000

/* Number of steps to take */
#define NUM_STEPS 300

int left_or_right(void)
{
	int r = rand();
	if (r < RAND_MAX/2) return -1;
	return +1;
}

int walk(int num_steps)
{
	int pt = 0;
	for (int i=0; i< num_steps; i++)
		pt += left_or_right();

	return pt;
}

int main(int argc, char* argv[])
{
	for (int ns=0; ns < NSAMP; ns++)
	{
		int pt = walk(NUM_STEPS);
		double x = pt / sqrt(NUM_STEPS);
	}
}
