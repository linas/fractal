/*
 * present.c
 *
 * hunt out group presentations
 *
 * Linas February 2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct _matrix {
	int dim;
	int *elts;
} Matrix;

Matrix * matrix_new (int dim)
{
	Matrix *mat = (Matrix *)malloc (sizeof (Matrix));
	mat->dim = dim;
	mat->elts = (int *)malloc (dim*dim*sizeof(int));
}

/* element acccessor */
#define MELT(mat,m,n) ((mat)->elts[(m)*(mat)->dim + (n)])

void matrix_set (Matrix *mat, int m, int n, int val)
{
	MELT(mat,m,n) = val;
}

void matrix_mult (Matrix *prod, Matrix *a, Matrix *b)
{
	int m,n,p;
	for (m=0; m<a->dim; m++)
	{
		for (n=0; n<a->dim; n++)
		{
			int acc = 0;
			for (p=0; p<a->dim; p++)
			{
				acc += MELT(a,m,p) * MELT(b, p,n);
			}
			MELT(prod,m,n) = acc;
		}
	}
}

int matrix_equal (Matrix *a, Matrix *b)
{
	int m,n;
	for (m=0; m<a->dim; m++)
	{
		for (n=0; n<a->dim; n++)
		{
			if (MELT(a,m,n) != MELT(b, m,n)) return 0;
		}
	}
	return 1;
}

/* ---------------------------------------------------- */
/* linked list of matrixes */
typedef struct _mlist MatList;

struct _mlist {
	MatList *next;
	Matrix *mat;
};

MatList *matlist_prepend (MatList *ml, Matrix *mat)
{
	MatList *node = (MatList *)malloc(sizeof (MatList));
	node->mat = mat;
	node->next = ml;

	return node;
}

/* Find a matching matrix in the list */
MatList *matlist_find (MatList *ptr, Matrix *mat)
{
	while (ptr)
	{
		if (matrix_equal (ptr->mat, mat)) return ptr;
		ptr = ptr->next;
	}
}

main ()
{
}
