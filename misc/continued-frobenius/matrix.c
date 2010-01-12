
/*
 * matrix.c
 *
 * Matrix utils
 *
 * Linas Jan 2004
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void
multmatrix (matrix *prod, matrix *ml, matrix *mr)
{
	int i,j,k;
	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			long double acc = 0.0;
			for (k=0; k<MS; k++)
			{
				acc += (*ml)[i][k] * (*mr)[k][j];
			}
			(*prod)[i][j] = acc;
		}
	}
}

void
identmatrix (matrix *e, long double val)
{
	int i,j;
	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			(*e)[i][j] = 0.0L;
		}
		(*e)[i][i] = val;
	}
}

void
copymatrix (matrix *to, matrix *from)
{
	int i,j;

	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			(*to)[i][j] = (*from)[i][j];
		}
	}
}

void
addmatrix (matrix *to, matrix *afrom, matrix *bfrom)
{
	int i,j;

	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			(*to)[i][j] = (*afrom)[i][j] +(*bfrom)[i][j];
		}
	}
}

void
lammatrix (matrix *to, long double lambda)
{
	int i;

	for (i=0; i<MS; i++)
	{
		(*to)[i][i] -= lambda;
	}
}

void
scalematrix (matrix *to, long double scale)
{
	int i,j;

	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			(*to)[i][j] *= scale;
		}
	}
}

