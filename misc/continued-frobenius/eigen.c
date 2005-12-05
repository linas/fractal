
/*
 * eigen.c
 *
 * Find the eigenvalues (zeros of the characteriztic polynomial i.e.
 * zeros of the determinant) of H.  H is the taylor-expanded matrix for 
 * the continued fraction map, expanded at y=1.
 *
 * Linas Jan 2004
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ache.h"
#include "zetafn.h"

// =================================================================

#define MS 100 // matrix dimension

typedef long double matrix[MS][MS];
typedef long double vector[MS];

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

/* This is a quick-n-dirty implementation for computing the determinant
 * of a matrix.
 */

typedef long double dubba;

/* alternate() returns the determinant of a submatrix of matrix m.
 *   vlen: input, size of submatrix within m
 *   ir, ic: vectors containing row numbers and column numbers
 *         from which the alternating product should be composed.
 *
 *  The implementation is recursive: the alternating product for
 *  larger matrices is expressed in terms of smaller ones.  The
 *  cases for the first 5 orders are computed explicitly, to minimize
 *  the overhead of a recursive call for these smaller sizes.  
 *
 *  The code for the first 5 cases was generated automatically
 *  to ensure correctness, i.e. to avoid typing mistakes.  See
 *  the ifdef'd portion below for the generator.
 *  
 *  Note that execution time grows like N factorial. Thus n=14
 *  might take less than an hour, n=15 might take more, and n=16 
 *  might take days.
 */
dubba
alternate (matrix *m, vector *ir, vector *ic, int vlen)
{
	dubba det;
	if (2 == vlen)
	{
		int ra, rb, ca, cb;
		ra = (*ir)[0];
		rb = (*ir)[1];

		ca = (*ic)[0];
		cb = (*ic)[1];

		det = (*m)[ra][ca] * (*m)[rb][cb];
		det -= (*m)[ra][cb] * (*m)[rb][ca];
		return det;
	}
	else if (3 == vlen)
	{
		int r0, r1, r2, c0, c1, c2;
		r0 = (*ir)[0];
		r1 = (*ir)[1];
		r2 = (*ir)[2];

		c0 = (*ic)[0];
		c1 = (*ic)[1];
		c2 = (*ic)[2];

		det = 
		(  ((*m)[r0][c0]) * (((*m)[r1][c1]) * ((*m)[r2][c2]) - ((*m)[r1][c2]) * ((*m)[r2][c1]))
		 - ((*m)[r0][c1]) * (((*m)[r1][c0]) * ((*m)[r2][c2]) - ((*m)[r1][c2]) * ((*m)[r2][c0]))
		 + ((*m)[r0][c2]) * (((*m)[r1][c0]) * ((*m)[r2][c1]) - ((*m)[r1][c1]) * ((*m)[r2][c0]))
		);
		return det;
	}
	else if (4 == vlen)
	{
		int r0, r1, r2, r3, c0, c1, c2, c3;
		r0 = (*ir)[0];
		r1 = (*ir)[1];
		r2 = (*ir)[2];
		r3 = (*ir)[3];

		c0 = (*ic)[0];
		c1 = (*ic)[1];
		c2 = (*ic)[2];
		c3 = (*ic)[3];

		det = 
		(  ((*m)[r0][c0]) * (  ((*m)[r1][c1]) * (((*m)[r2][c2]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c2]))
		 - ((*m)[r1][c2]) * (((*m)[r2][c1]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c1]))
		 + ((*m)[r1][c3]) * (((*m)[r2][c1]) * ((*m)[r3][c2]) - ((*m)[r2][c2]) * ((*m)[r3][c1]))
		)
		 - ((*m)[r0][c1]) * (  ((*m)[r1][c0]) * (((*m)[r2][c2]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c2]))
		 - ((*m)[r1][c2]) * (((*m)[r2][c0]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c0]))
		 + ((*m)[r1][c3]) * (((*m)[r2][c0]) * ((*m)[r3][c2]) - ((*m)[r2][c2]) * ((*m)[r3][c0]))
		)
		 + ((*m)[r0][c2]) * (  ((*m)[r1][c0]) * (((*m)[r2][c1]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c1]))
		 - ((*m)[r1][c1]) * (((*m)[r2][c0]) * ((*m)[r3][c3]) - ((*m)[r2][c3]) * ((*m)[r3][c0]))
		 + ((*m)[r1][c3]) * (((*m)[r2][c0]) * ((*m)[r3][c1]) - ((*m)[r2][c1]) * ((*m)[r3][c0]))
		)
		 - ((*m)[r0][c3]) * (  ((*m)[r1][c0]) * (((*m)[r2][c1]) * ((*m)[r3][c2]) - ((*m)[r2][c2]) * ((*m)[r3][c1]))
		 - ((*m)[r1][c1]) * (((*m)[r2][c0]) * ((*m)[r3][c2]) - ((*m)[r2][c2]) * ((*m)[r3][c0]))
		 + ((*m)[r1][c2]) * (((*m)[r2][c0]) * ((*m)[r3][c1]) - ((*m)[r2][c1]) * ((*m)[r3][c0]))
		));
		return det;
	}
	else if (5 == vlen)
	{
		int r0, r1, r2, r3, r4, c0, c1, c2, c3, c4;
		r0 = (*ir)[0];
		r1 = (*ir)[1];
		r2 = (*ir)[2];
		r3 = (*ir)[3];
		r4 = (*ir)[4];

		c0 = (*ic)[0];
		c1 = (*ic)[1];
		c2 = (*ic)[2];
		c3 = (*ic)[3];
		c4 = (*ic)[4];

		det = 
		(  ((*m)[r0][c0]) * (  ((*m)[r1][c1]) * (  ((*m)[r2][c2]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		)
		 - ((*m)[r1][c2]) * (  ((*m)[r2][c1]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		)
		 + ((*m)[r1][c3]) * (  ((*m)[r2][c1]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		)
		 - ((*m)[r1][c4]) * (  ((*m)[r2][c1]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		)
		)
		 - ((*m)[r0][c1]) * (  ((*m)[r1][c0]) * (  ((*m)[r2][c2]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		)
		 - ((*m)[r1][c2]) * (  ((*m)[r2][c0]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		)
		 + ((*m)[r1][c3]) * (  ((*m)[r2][c0]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		)
		 - ((*m)[r1][c4]) * (  ((*m)[r2][c0]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		)
		)
		 + ((*m)[r0][c2]) * (  ((*m)[r1][c0]) * (  ((*m)[r2][c1]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		)
		 - ((*m)[r1][c1]) * (  ((*m)[r2][c0]) * (((*m)[r3][c3]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c3]))
		 - ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		)
		 + ((*m)[r1][c3]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)
		 - ((*m)[r1][c4]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)
		)
		 - ((*m)[r0][c3]) * (  ((*m)[r1][c0]) * (  ((*m)[r2][c1]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		)
		 - ((*m)[r1][c1]) * (  ((*m)[r2][c0]) * (((*m)[r3][c2]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		)
		 + ((*m)[r1][c2]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c4]) - ((*m)[r3][c4]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c4]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)
		 - ((*m)[r1][c4]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)
		)
		 + ((*m)[r0][c4]) * (  ((*m)[r1][c0]) * (  ((*m)[r2][c1]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		)
		 - ((*m)[r1][c1]) * (  ((*m)[r2][c0]) * (((*m)[r3][c2]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c2]))
		 - ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		)
		 + ((*m)[r1][c2]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c3]) - ((*m)[r3][c3]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c3]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)
		 - ((*m)[r1][c3]) * (  ((*m)[r2][c0]) * (((*m)[r3][c1]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c1]))
		 - ((*m)[r2][c1]) * (((*m)[r3][c0]) * ((*m)[r4][c2]) - ((*m)[r3][c2]) * ((*m)[r4][c0]))
		 + ((*m)[r2][c2]) * (((*m)[r3][c0]) * ((*m)[r4][c1]) - ((*m)[r3][c1]) * ((*m)[r4][c0]))
		)));
		return det;
	}
	else
	{
		vector row, col;

		int rrr = (*ir)[0];
		
		// skip the top row when working submatrices
		int i,j;
		for (i=0; i<vlen-1; i++)
		{
			row[i] = (*ir)[i+1];
		}
		
		dubba sign = 1.0L;
		dubba acc = 0.0L;
		for (i=0; i<vlen; i++)
		{
			int k=0;
			for (j=0; j<vlen-1; j++)
			{
				if (k==i) k++;  // skip this column
				col[j] = (*ic)[k];
				k++;
			}

			int ccc = (*ic)[i];
			acc += sign * (*m)[rrr][ccc] * alternate (m, &row, &col, vlen-1);
			sign = -sign;
		}
		return acc;
	}
}

/* Return determinant of matrix m, of dimension dim */

dubba
determinant (matrix *m, int dim)
{
	int i;
	vector row, col;

	for (i=0; i<dim; i++)
	{
		row[i] = i;
		col[i] = i;
	}

	return alternate (m, &row, &col, dim);
}

#ifdef TEST
main ()
{
	matrix m;
	int  i, j;
	for (i=0; i<5; i++) 
	{
		for (j=0; j<5; j++)
		{
			// m[i][j] = j*i;
			m[i][j] = 0.0;
		}
		m[i][i] = 1.0;
	}

	printf ("its %Lg\n", determinant (&m, 5));
	return 0;
}
#endif

/* ===================================================================== */
#ifdef DET_GENERATOR

typedef int vector[30];

/* This code generates the source code for the lower orders.
 * Its straightforward to audit
 */
void
alt (vector *ir, vector *ic, int vlen)
{
	if (2 == vlen)
	{
		int ra, rb, ca, cb;
		ra = (*ir)[0];
		rb = (*ir)[1];

		ca = (*ic)[0];
		cb = (*ic)[1];

		printf ("(((*m)[r%d][c%d]) * ((*m)[r%d][c%d]) - ((*m)[r%d][c%d]) * ((*m)[r%d][c%d]))\n",
			ra,ca,rb,cb, ra,cb,rb,ca);
	}
	else
	{
		vector row, col;

		int rrr = (*ir)[0];
		int i,j;
		for (i=0; i<vlen-1; i++)
		{
			row[i] = (*ir)[i+1];
		}

		int sign = 1;
		printf ("(  ");
		for (i=0; i<vlen; i++)
		{
			int k=0;
			for (j=0; j<vlen-1; j++)
			{
				if (k==i) k++;  // skip this column
				col[j] = (*ic)[k];
				k++;
			}

			int ccc = (*ic)[i];
			if (i != 0) {
				if (sign > 0) printf (" + ");
				else printf (" - ");
			}
			printf ("((*m)[r%d][c%d]) * ", rrr, ccc);
			alt (&row, &col, vlen-1);
			sign = -sign;
		}
		printf (")\n");
	}
}

void
det_gen (int dim)
{
	int i;
	vector row, col;

	for (i=0; i<dim; i++)
	{
		row[i] = i;
		col[i] = i;
	}

	alt (&row, &col, dim);
}

main ()
{
	det_gen (5);
	return 0;
}
#endif

/* ============================================================= */

static matrix h;

void
setup(void)
{
	int i,j;
	for (i=0; i<40; i++)
	{
		for (j=0; j<40; j++)
		{
			h[i][j] = ache_mp(i,j);
		}
	}
}

long double get_det (long double lambda, int dim)
{
	matrix work;
	copymatrix (&work, &h);
	lammatrix (&work, lambda);

	long double det = determinant (&work, dim);
	return det;
}

/* ============================================================= */

#define UPDATE   \
		if (0.0 < fa*fc)                 \
		{                 \
			have_prev_a = 1;                \
			prev_a = a;                \
			prev_fa = fa;                \
			a=c; fa = fc;                 \
		}                 \
		else                 \
		{                 \
			have_prev_b = 1;                \
			prev_b = b;                \
			prev_fb = fb;                \
			b=c; fb=fc;                 \
		}


#define TEST_DONE(STR)  \
	printf ("# " STR " its a=%15.12Lg b=%15.12Lg fc=%Lg delt=%Lg\n", a,b, fc, b-a); \
	printf ("# date: %s", ctime (({time_t now=time(0); &now;})));  \
	if (fmax > fabsl(fc)) break;   \
	if (1.0e-8 > fabsl(a-b)) break;   \
	fflush (stdout);  


long double 
zero (int dim, long double a, long double b)
{
	long double fa, fb, c, fc;
	int have_prev_a = 0;
	int have_prev_b = 0;
	long double prev_a=0.0, prev_b=0.0, prev_fa=0.0, prev_fb=0.0;

	if (a>b)
	{
		long double tmp = a;
		a=b;
		b=tmp;
	}
	
	fa = get_det (a, dim);
	printf ("# a=%Lg fa=%15.12Lg\n", a, fa);
	printf ("# date: %s", ctime (({time_t now=time(0); &now;})));
	fflush (stdout);
	fb = get_det (b, dim);
	printf ("# b=%Lg fb=%15.12Lg\n", b, fb);
	printf ("# date: %s", ctime (({time_t now=time(0); &now;})));
	fflush (stdout);

	long double fmax = 1.0e-5 * fabs (fb-fa);
	printf ("# look for f-bound=%Lg\n", fmax);

	while (1)
	{
		// linear interpolate
		c = (-fa*b + fb*a) / (fb-fa);
		fc = get_det(c, dim);
		UPDATE;
		TEST_DONE("interp")

		// attempt to linear extrapolate ... else bisect
		if (have_prev_a)
		{
			c = (-fa*prev_a + prev_fa*a) / (prev_fa-fa);
			if ((a<c) && (c<b))
			{
				fc = get_det(c, dim);
				UPDATE;
				TEST_DONE("extra-a")
				continue;
			}
		}
		else
		if (have_prev_b)
		{
			c = (-fb*prev_b + prev_fb*b) / (prev_fb-fb);
			if ((a<c) && (c<b))
			{
				fc = get_det(c, dim);
				UPDATE;
				TEST_DONE("extra-b")
				continue;
			}
		}
		
		// if we got here, extrapolation failed, so bisect
		c = 0.5*(a+b);
		fc = get_det(c, dim);
		UPDATE;
		TEST_DONE("bisect")
	}
	
	// for the final answer, linear interpolate
	c = (-fa*b + fb*a) / (fb-fa);
	return c;
}


/* ============================================================= */

#define DO_ZERO
#ifdef DO_ZERO
int main ()
{
	setup();

	printf ("# \n");
	printf ("# Eigenvalues of truncated H matrix\n");
	printf ("# Solves characteristic equation for zeros\n");
	printf ("# output: location of zero as rank of matrix\n");
	printf ("# date: %s", ctime (({time_t now=time(0); &now;})));
	int dim;
	for (dim=16; dim<30; dim++)
	{
		printf ("# begin zero finding for dim=%d\n", dim);
		fflush (stdout);
		// to fit lam-1, start at dim=8 and bound as below
		// long double zz = zero (dim, -0.29, -0.31);
		// 
		// to fit lam-1, start at dim=10 and bound as below
		// long double zz = zero (dim, -0.303, -0.304);
		// 
		// to fit lam-1, start at dim=16 and bound as below
		// long double zz = zero (dim, -0.303651584, -0.30367);
		// 
		// to fit lam-1, start at dim=17 and bound as below
		// long double zz = zero (dim, -0.30365759, -0.303662);
		// 
		// to fit lam-2, start at dim 8 and bound as below
		// long double zz = zero (dim, 0.092, 0.105);
		// 
		// to fit lam-2, start at dim 10 and bound as below
		// long double zz = zero (dim, 0.0995, 0.102);
		// 
		// to fit lam-2, start at dim 16 and bound as below
		// warning this will take weeks of cpu time
		long double zz = zero (dim, 0.10085, 0.1009);
		// 
		// to fit lam-3, start at dim=8 and bound as below
		// long double zz = zero (dim, -0.028, -0.037);
		//
		// to fit lam-3, start at dim=10 and bound as below
		// long double zz = zero (dim, -0.033, -0.0363);
		//
		// to fit lam-4, start at dim=10 and bound as below
		// long double zz = zero (dim, 0.010, 0.013);
		// 
		// to fit lam-5, start at dim=10 and bound as below
		// long double zz = zero (dim, -0.0019, -0.005);
		//
		// to fit lam-5, start at dim=15 and bound as below
		// long double zz = zero (dim, -0.00436, -0.0046);
		//
		// to fit lam-6, start at dim=10 and bound as below
		// long double zz = zero (dim, 0.0002, 0.0017);
		//
		// to fit lam-6, start at dim=13 and bound as below
		// long double zz = zero (dim, 0.00084, 0.0017);
		// 
		// to fit lam-7, start at dim=11 and bound as below
		//long double zz = zero (dim, -0.00003, -0.0005);
		//
		// to fit lam-7, start at dim=12 and bound as below
		// long double zz = zero (dim, -0.0001, -0.0005);
		// 
		// to fit lam-7, start at dim=14 and bound as below
		// long double zz = zero (dim, -0.0002, -0.0004);
		//

		printf ("%d\t%15.12Lg\n", dim, zz);
		fflush (stdout);
	}
	return 0;
}
#endif

//  #define DO_GRAPHS
#ifdef DO_GRAPHS
int
main(int argc, char *argv[])
{
	setup();

	int dim = 3;
	if (argc == 2) dim = atoi(argv[1]);
	printf ("# \n");
	printf ("# characterisitic equation \n");
	printf ("# Dimension=%d\n", dim);

	long double lambda;
	for (lambda = 0.05; lambda > -0.05; lambda -= 0.0005)
	{
		//long double det = get_det (lambda, dim);
		// printf ("its %Lg    %Lg\n", lambda, det);
#define MULTI
#ifdef MULTI
		printf ("%Lg", lambda);
		for (dim=2; dim<11; dim++)
		{
			long double det = get_det (lambda, dim);
			printf ("\t%Lg", det);
		}
		printf ("\n");

#endif
		fflush (stdout);
	}
	
	return 0;
}
#endif
