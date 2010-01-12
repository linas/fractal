
/*
 * matrix.c
 *
 * Matrix utils
 *
 * Linas Jan 2004
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>o
#include "matrix.h"

void
multmatrix (int dim, matrix *prod, matrix *ml, matrix *mr)
{
	int i,j,k;
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			long double acc = 0.0;
			for (k=0; k<dim; k++)
			{
				acc += (*ml)[i][k] * (*mr)[k][j];
			}
			(*prod)[i][j] = acc;
		}
	}
}

void
identmatrix (int dim, matrix *e, long double val)
{
	int i,j;
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			(*e)[i][j] = 0.0L;
		}
		(*e)[i][i] = val;
	}
}

void
copymatrix (int dim, matrix *to, matrix *from)
{
	int i,j;

	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			(*to)[i][j] = (*from)[i][j];
		}
	}
}

void
addmatrix (int dim, matrix *to, matrix *afrom, matrix *bfrom)
{
	int i,j;

	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			(*to)[i][j] = (*afrom)[i][j] +(*bfrom)[i][j];
		}
	}
}

void
lammatrix (int dim, matrix *to, long double lambda)
{
	int i;

	for (i=0; i<dim; i++)
	{
		(*to)[i][i] -= lambda;
	}
}

void
scalematrix (int dim, matrix *to, long double scale)
{
	int i,j;

	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			(*to)[i][j] *= scale;
		}
	}
}

