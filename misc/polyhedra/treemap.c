
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

typedef struct {
	int cnt;
	double sum;
	double sumsq;
} Stats;


typedef struct _vertex Vertex;
typedef struct _edge Edge;

/* Directed edge, from from to to */
struct _edge {
	Vertex *from;
	Vertex *to;
	Edge *left;
	Edge *right;
	Stats stats;
};

struct _vertex {
	int id;
	Edge edge[3];  // outgoing, directed edge
	Stats stats;
};

Vertex * tetrahedron_new (void)
{
	int i, j;
	Vertex * t = (Vertex *) malloc (4*sizeof (Vertex));

	for (i=0; i<4; i++)
	{
		t[i].id = i;
		t[i].stats.cnt = 0;
		t[i].stats.sum = 0.0;
		t[i].stats.sumsq = 0.0;

		for (j=0; j<2; j++)
		{
			t[i].edge[j].stats.cnt = 0;
			t[i].edge[j].stats.sum = 0.0;
			t[i].edge[j].stats.sumsq = 0.0;
			t[i].edge[0].to = &t[i];
			t[i].edge[1].from = &t[i];
			t[i].edge[2].from = &t[i];
		}
	}

	t[0].edge[0].left  = &t[1].edge[0];
	t[0].edge[0].right = &t[2].edge[0];
	t[0].edge[0].from  = &t[3];

	t[0].edge[1].to    = &t[1];
	t[0].edge[1].left  = &t[1].edge[1];
	t[0].edge[1].right = &t[1].edge[2];

	t[0].edge[2].to    = &t[2];
	t[0].edge[2].left  = &t[2].edge[1];
	t[0].edge[2].right = &t[2].edge[2];

	/* hard */
	t[0].edge[0].from  = &t[3];
	t[0].edge[1].to    = &t[1];
	t[0].edge[2].to    = &t[2];

	for (i=0; i<4; i++)
	{
		t[i].edge[0].left  = &t[i].edge[1].to ->edge[0];  // t[0].edge[0].left  = &t[1].edge[0];
		t[i].edge[0].left  = &t[i].edge[2].to ->edge[0];  //  t[i].edge[0].right = &t[2].edge[0];
	
		t[i].edge[1].left  = &t[1].edge[1];
		t[i].edge[1].right = &t[1].edge[2];
	
		t[i].edge[2].left  = &t[2].edge[1];
		t[i].edge[2].right = &t[2].edge[2];
	}
	/*

	t[1].left  = &t[3];
	t[1].right = &t[2];

	t[2].left  = &t[3];
	t[2].right = &t[0];

	t[3].left  = &t[1];
	t[3].right = &t[0];
*/

}

main ()
{
	tetrahedron_new();
}
