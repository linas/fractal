/*
 * random-walk.c
 *
 * Gaussian random walk.
 * Linas Vepstas 24 March 2025
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

/* Number of histogram bins. */
#define HISTO_NBINS 302

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

	double norm = ((double) (HISTO_NBINS-2)) / ((double) tot_cnt);

	for (int i=0; i< HISTO_NBINS; i++)
	{
		double x = (i-1) / ((double) (HISTO_NBINS-2));
		double y = norm * histogram[i];
		printf("%d	%f	%d	%f\n", i, x, histogram[i], y);
	}
	fflush(stdout);
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s nsteps nsamp\n", argv[0]);
		exit(-1);
	}

	/* Number of samples */
	int nsamp = atoi(argv[2]);

	/* Number of steps to take */
	int num_steps = atoi(argv[1]);
	double norm = 1.0 / sqrt(num_steps);

	/* Also missing sqrt 2pi but whatever. */

	/* six sigma */
	norm /= 6.0;

	histo_init();
	for (int ns=0; ns < nsamp; ns++)
	{
		int pt = walk(num_steps);
		double x = norm * pt;
		x += 0.5;
		histo_accum(x);
	}

	printf("#\n# num-steps=%d nsamp=%d\n#\n", num_steps, nsamp);
	histo_print();
}
