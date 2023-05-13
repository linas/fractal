
/*
 * fliptree.C
 *
 * show tree for 1 to infinitiy
 *
 * Linas November 2004
 */

#include <stdio.h>
#include <stdlib.h>

double * flip_tree (int depth)
{
	double * tree;
	
	int size = 2<<(depth-1);

	tree = (double *) malloc (size * sizeof (double));

	int pos = size/2 -1;

	printf ("# tree size=%d depth=%d\n", size-1, depth);
	tree[size-1] = 1.0;
	// printf ("pos for 2 is %d\n", pos);
	tree[pos] = 2.0;
	int lev;
	int tn = 2;
	for (lev=2; lev<=depth; lev++)
	{
		int delt = size/tn;
		int step = delt/2;
		int start = step -1;
		int end = size-1 - step;
		tree[start] = 2*tn;
		// printf ("lve=%d   lowb=%d val=%g\n", lev, start, tree[start]);
		tree[end] =  0.5 + 0.5*tree[end-step];
		// printf ("lev=%d hibb=%d val=%g\n", lev, end, tree[end]);

		// printf ("start at %d walk by %d end at %d\n", start, delt, end);
		int i;
		for (i=start+delt; i<end; i+= delt)
		{
			tree[i] = 0.5 * (tree[i-step] + tree[i+step]);
			// printf ("do i=%d, got val=%g\n", i, tree[i]);
		}
		tn *= 2;
	}

	return tree;
}

main ()
{
	double * tree;

	int level = 6;
	tree = flip_tree (level);

	int sz = 1<<level;
	int i;
	for (i=0; i<sz-1; i++)
	{
		double x = (i+1)/((double) sz);
		double slope = tree[i+1]-tree[i];
		slope *= (double) sz;
		printf ("%d %g  \t1/x=%g  \t%g    \tslope=%g\n", i, x, 1.0/x, tree[i], slope);
	}
}
