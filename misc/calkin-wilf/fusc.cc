/*
 * Calkin-Wilf tree fusc function
 *
 * September 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "question.h"

/**
 * The "fusc" fuction from the Calkin-Wilf tree.
 * This gives the denominator of the stern-brocot tree.
 *
 * There is nothing (strictly) analogous to it for the numerator
 * of the stern-brocot tree, because of problematic values at
 * exact powers-of-two, which must be either 0 or 1 depending
 * on where in the sequence we are.
 */
unsigned long fusc(unsigned long n)
{
	if (0 == n) return 0;
	if (1 == n) return 1;
	unsigned long m = n >> 1;
	if (0 == n%2) return fusc(m);
	return fusc(m) + fusc(m+1);
}

unsigned long nusc(unsigned long n)
{
	if (0 == n) return 0;
	if (1 == n) return 0;
	if (2 == n) return 1;
	unsigned long m = n >> 1;
	if (0 == n%2) return nusc(m);
	return nusc(m+1) + nusc(m);
}

unsigned long section_sum(int section)
{
	unsigned long bot = 1<<section;
	unsigned long top = 2*bot;

	unsigned long sum = 0;
	for (unsigned long i=bot; i<= top; i++)
	{
		unsigned long fu = fusc(i);
		sum += fu;
	}
	return sum;
}

void print_bro(int section)
{
	unsigned long bot = 1<<section;
	unsigned long top = 2*bot;

	unsigned long p,q;
	for (unsigned long i=bot; i<= top; i++)
	{
		unsigned long fu = fusc(i);
		unsigned long nu = nusc(i);
		stern_brocot_tree(i-bot, section, &p, &q);
		printf("%lu	%lu	%lu	%lu	%lu\n", i-bot, fu, nu, p, q);
		if (0 != fu-q)
		{
			printf("oooooooo nooooo!\n");
			exit(1);
		}
	}
}

void print_sum(int section)
{
	unsigned long bot = 1<<section;
	unsigned long top = 2*bot;

	unsigned long sum = 0;
	double norm = 1.0 / (double) section_sum(section);

	for (unsigned long i=bot; i<= top; i++)
	{
		unsigned long fu = fusc(i);
		double x = (i - bot) / ((double) bot);
		sum += fu;
		double rsum = sum * norm;
		printf("%lu	%g	%lu	%lu	%g\n", i, x, fu, sum, rsum);
	}
}

void print_weighted_sum(int section)
{
	unsigned long bot = 1<<section;
	unsigned long top = 2*bot;

	unsigned long sum = 0;
	double wsum = 0;
	double norm = 1.0 / (double) section_sum(section);

	unsigned long p, q;
	double prev = 0.0;

	for (unsigned long i=bot+1; i<= top; i++)
	{
		unsigned long fu = fusc(i);
		double x = (i - bot) / ((double) bot);
		sum += fu;
		double rsum = sum * norm;

		stern_brocot_tree(i-bot, section, &p, &q);
		double qinv = ((double) p)/ ((double) q);
		double delta = qinv - prev;
		prev = qinv;
		wsum += fu / delta;
		printf("%lu	%g	%g	%g	%lu	%lu	%g	%g\n", i-bot, x, qinv, delta, fu, sum, rsum, wsum);
	}
}

int main(int argc, char * argv[])
{
	int section = atof(argv[1]);

#if 0
	unsigned long ssum = 0;
	for (unsigned long i=0; i<= 30; i++)
	{
		unsigned long s = section_sum(i);
		ssum += s;
		printf("its %lu %lu	%lu\n",i, s, ssum);
	}
#endif

	// print_sum(section);
	// print_weighted_sum(section);
	print_bro(section);
}
