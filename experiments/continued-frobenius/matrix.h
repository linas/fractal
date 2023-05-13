
/*
 * matrix.h
 *
 * Matrix utils
 *
 * Linas Jan 2004
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef long double matrix[MS][MS];
typedef long double vector[MS];

static inline void
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

static inline void
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

static inline void
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

static inline void
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

static inline void
lammatrix (matrix *to, long double lambda)
{
	int i;

	for (i=0; i<MS; i++)
	{
		(*to)[i][i] -= lambda;
	}
}

static inline void
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


#endif /* __MATRIX_H__ */
