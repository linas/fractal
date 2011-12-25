
/*
 * julia.c
 * simpel julia set explorer for the mandelbrot bulb
 *
 * Linas Vepstas December 2011
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void julia_tree (complex double * tree, int depth, complex double c, complex double w)
{
	complex double z = w - c;
	double r = sqrt(cabs (z));
	double theta = 0.5 * carg(z);

	complex double v = r * cexp (I*theta);

	int shift = 1<<(depth-1);
   if (1 == depth)
	{
   	tree[0] = v;
   	tree[1] = -v;
		return;
	}

// printf("duude v %g and %g\n", creal(v), cimag(v));
   julia_tree(tree, depth-1, c, v);
   julia_tree(&tree[shift], depth-1, c, -v);

}

complex double * make_tree (int depth, complex double c)
{
   // depth is same as width in bits.
	int tree_width = 1 << depth;
	complex double *tree = malloc (tree_width * sizeof (complex double));

	julia_tree(tree, depth, c, 0.0);

	return tree;
}

int main (int argc, char * argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <tree-depth>\n", argv[0]);
		exit(1);
	}

	int depth = atoi (argv[1]);

	complex double * tree;
	tree = make_tree(depth, -0.8+0.6*I);

	int sz = 1<< depth;
	int i;
	for (i=0; i<sz; i++)
	{
		complex double v = tree[i];
		printf("%d\t%g\t%g\n", i, creal(v), cimag(v));
	}
}
