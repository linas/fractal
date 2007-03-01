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
#include <string.h>

typedef struct _matrix {
	int dim;
	int *elts;
} Matrix;

Matrix * matrix_new (int dim)
{
	Matrix *mat = (Matrix *)malloc (sizeof (Matrix));
	mat->dim = dim;
	mat->elts = (int *)malloc (dim*dim*sizeof(int));
	return mat;
}

/* element acccessor */
#define MELT(mat,m,n) ((mat)->elts[(m)*(mat)->dim + (n)])

void matrix_set (Matrix *mat, int m, int n, int val)
{
	MELT(mat,m,n) = val;
}

void matrix_copy (Matrix *to, Matrix *from)
{
	int m,n;
	for (m=0; m<from->dim; m++)
	{
		for (n=0; n<from->dim; n++)
		{
			MELT(to,m,n) = MELT(from, m,n);
		}
	}
}

Matrix * matrix_dup (Matrix *from)
{
	Matrix *mat = matrix_new (from->dim);
	matrix_copy (mat, from);
	return mat;
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
	char *str;
};

MatList *matlist_prepend (MatList *ml, Matrix *mat, char *str, char letter)
{
	MatList *node = (MatList *)malloc(sizeof (MatList));
	node->mat = mat;
	node->next = ml;
	int len = strlen(str);
	node->str = malloc (len+2);
	strcpy(node->str, str);
	node->str[len] = letter;
	node->str[len+1] = 0x0;

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
	return NULL;
}

/* ---------------------------------------------------- */
/* master struct */

typedef struct _present {
	MatList *generators;
	MatList *words;
} Present;

void walk_tree (Present *pr, MatList *node, int depth)
{
	if (0 == depth) return;
	depth --;

	MatList *gen = pr->generators;
	while (gen)
	{
		Matrix *mat = matrix_new (node->mat->dim);
		matrix_mult (mat, node->mat, gen->mat);

		pr->words = matlist_prepend (pr->words, mat, node->str, gen->str[0]);	
		gen = gen->next;
	}
}

/* ---------------------------------------------------- */

Present * setup_b3 (void)
{
	Present *pr = (Present *) malloc (sizeof (Present));
	pr->generators = NULL;
	pr->words = NULL;


	/* identity matrix */
	Matrix *e = matrix_new (2);
	MELT(e, 0,0) = 1;
	MELT(e, 0,1) = 0;
	MELT(e, 1,0) = 0;
	MELT(e, 1,1) = 1;
	pr->words = matlist_prepend (NULL, e, "", 'E');
	
	MatList *ml = NULL;

	/* Set up sigma 1 and sigma 2 of braid group B_3 */
	Matrix *s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = 1;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'A');
	
	Matrix *s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = -1;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'B');
	
	/* and now thier inverses (by hand) */
	s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = -1;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'a');
	
	s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = 1;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'b');

	pr->generators = ml;
	return pr;
}

main ()
{
	Present *pr = setup_b3();
	walk_tree (pr, pr->words, 3);
}
