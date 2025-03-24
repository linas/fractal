/*
 * random-walk.c
 *
 * Gaussian random walk.
 * Linas Vepstas 24 March 2025
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Number of histogram bins. */
#define HISTO_NBINS 302

/* Number of samples */
#define NSAMP 5000

/* Number of steps to take */
#define NUM_STEPS 3000

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

int histogram[HISTO_NBINS];

void histo_init(void)
{
	for (int i=0; i< HISTO_NBINS; i++)
		histogram[i] = 0;
}

void histo_accum(double x)
{
	int bin = 1 + floor(x * (HISTO_NBINS-2));
	if (bin < 0) bin = 0;
	if (HISTO_NBINS <= bin) bin = HISTO_NBINS-1;
	histogram[bin] += 1;
}

void histo_print(void)
{
	int tot_cnt = 0;
	for (int i=0; i< HISTO_NBINS; i++)
		tot_cnt += histogram[i];

	for (int i=0; i< HISTO_NBINS; i++)
	{
		double x = (i-1) / ((double) (HISTO_NBINS-2));
		double y = histogram[i] / ((double) tot_cnt);
		printf("%d	%f	%d	%f\n", i, x, histogram[i], y);
	}
	fflush(stdout);
}

int main(int argc, char* argv[])
{
	histo_init();
	for (int ns=0; ns < NSAMP; ns++)
	{
		int pt = walk(NUM_STEPS);
		// double x = pt / sqrt(NUM_STEPS);
		double x = pt / ((double) NUM_STEPS);
		x += 0.5;
		histo_accum(x);
	}

	histo_print();
}
