
/*
 * lr-tree.c
 *
 * left-right tree. A finite length string of L's and R's encodes
 * a series of left-right moves on the binary tree. Idea is to
 * label nodes on the tree with dyadics: so, no-move=1/2 and
 * L=1/4 and R=3/4 and so on.
 *
 * String of LR's are encoded in bits of real number 0 <= x <= 1
 * justified to decimal point (i.e. as dyadic).
 */

#include <stdio.h>
#include <stdlib.h>

/* tree_to_dyadic -- return dydic value of tree moves 
 * @x: real-encoded string
 * @len: string len
 */

double tree_to_dyadic (double x, int len)
{
	double label=0.5;
	int i;

	double step=0.25;
	for (i=0; i<len; i++)
	{
		if (0.5 <= x)
		{
			label += step;
			x = 2.0*x -1.0;
		}
		else
		{
			label -= step;
			x = 2.0*x;
		}

		step *= 0.5;
	}

	return label;
}

main ()
{
	int k,n, len;

	n=10;

	for (len=0; len<n; len++)
	{
		int m = 1<<len;
		for (k=0; k<m; k++)
		{
			// double x = ((double) k)/ ((double) m);
			double x = ((double) k)/ ((double) m);
			double y = tree_to_dyadic (x,len);

			printf ("%d	%d	%8.6g	%8.6g\n",k,len,x,y);
		}
	}
}
