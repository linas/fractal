
/*
 * treemap.c
 *
 * Map a binary tree onto a polyhedron
 *
 * Linas February 2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct _vertex Vertex;

struct _vertex {
	Vertex *left;
	Vertex *right;
	int cnt;
	double sum;
	double sumsq;
};

Vertex * tetrahedron_new (void)
{
	Vertex * t = (Vertex *) malloc (4*sizeof (Vertex));

	t[0].left = &t[1];
	t[0].right = &t[2];

}

main ()
{
	tetrahedron_new();
}
