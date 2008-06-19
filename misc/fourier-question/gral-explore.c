/*
 * Explore matrix elements of fourier of minkowski question mark 
 *
 * Linas June 2008
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gral-simple.h"

void slice_row(int row, int colmin, int colmax)
{
	int j;

	for (j=colmin; j<colmax; j++)
	{
		long double re = 0.0L;
		long double im = 0.0L;

		make_elt(row, j, &re, &im);
		printf("%d	%d	%Lg	%Lg\n", row, j, re, im);
	}
}

void slice_col(int rowmin, int rowmax, int col)
{
	int i;

	for (i=rowmin; i<rowmax; i++)
	{
		long double re = 0.0L;
		long double im = 0.0L;

		make_elt(i, col, &re, &im);
		printf("%d	%d	%Lg	%Lg\n", i, col, re, im);
	}
}


main(int argc, char * argv[])
{
	set_npts(123123);

	int row = atoi(argv[1]);
	slice_row(row, -20, 360);
	//
	int col = row;
	// slice_col(-20, 160, col);
}
