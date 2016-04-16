/*
 * Generating functions for greatest prime factors.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>

#include <gpf.h>

int main()
{
	for (int n=1; n<100; n++)
	{
		int g = gpf(n);
		printf("duuude n=%d gpf=%d\n", n, g);
	}
}
