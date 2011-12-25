
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

void julia_tree (complex double * tree, int depth,
                 complex double c, complex double w,
                 double rot)
{
	complex double z = w - c;
	double r = sqrt(cabs (z));
	double theta = 0.5 * carg(z);
	if (theta > rot) theta -= M_PI;
	if (theta < rot-M_PI) theta += M_PI;

	complex double v = r * cexp (I*theta);

	int shift = 1<<(depth-1);
   if (1 == depth)
	{
   	tree[0] = v;
   	tree[1] = -v;
		return;
	}

// printf("duude v %g and %g\n", creal(v), cimag(v));
   julia_tree(tree, depth-1, c, v, rot);
   julia_tree(&tree[shift], depth-1, c, -v, rot);

}

complex double * make_tree (int depth, complex double c, double rot)
{
   // depth is same as width in bits.
	int tree_width = 1 << depth;
	complex double *tree = malloc (tree_width * sizeof (complex double));

	julia_tree(tree, depth, c, 0.0, rot);

	return tree;
}

unsigned int bit_reverse (int depth, unsigned int str)
{
   int i;
	unsigned int rev = 0;

	for (i=0; i<depth; i++)
	{
		rev <<= 1;
		if (str & 0x1) rev |= 0x1;
		str >>= 1;
	}

	return rev;
}

int main (int argc, char * argv[])
{
	if (argc != 5)
	{
		fprintf(stderr, "Usage: %s <tree-depth> <re> <im> <rot>\n", argv[0]);
		exit(1);
	}

	int depth = atoi (argv[1]);
	double re = atof (argv[2]);
	double im = atof (argv[3]);
	double rot = atof (argv[4]);

	complex double * tree;
	tree = make_tree(depth, re+im*I, rot);

	unsigned int sz = 1<< depth;
	unsigned int i;
	for (i=0; i<sz; i++)
	{
		unsigned int rev = bit_reverse(depth, i);
		complex double v = tree[rev];
		double x = ((double) i) / ((double) sz);
		printf("%d	%d	%g	%g	%g\n", i, rev, x, creal(v), cimag(v));
	}
	double x = 1.0;
	complex double v = tree[0];
	printf("%d	%d	%g	%g	%g\n", sz, 0, x, creal(v), cimag(v));

	return 0;
}
