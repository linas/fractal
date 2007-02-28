
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
	int i, j, k,m;
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

#ifdef EXAMPLE_SETUP
	t[0].edge[0].left  = &t[1].edge[0];
	t[0].edge[0].right = &t[2].edge[2];
	t[0].edge[0].from  = &t[3];

	t[0].edge[1].to    = &t[1];
	t[0].edge[1].left  = &t[1].edge[1];
	t[0].edge[1].right = &t[1].edge[2];

	t[0].edge[2].to    = &t[2];
	t[0].edge[2].left  = &t[2].edge[1];
	t[0].edge[2].right = &t[2].edge[2];
#endif

	/* hard-coded connectivity */
	t[0].edge[0].from  = &t[3];
	t[0].edge[1].to    = &t[1];
	t[0].edge[2].to    = &t[2];

	t[1].edge[0].from  = &t[0];
	t[1].edge[1].to    = &t[3];
	t[1].edge[2].to    = &t[2];

	t[2].edge[0].from  = &t[1];
	t[2].edge[1].to    = &t[3];
	t[2].edge[2].to    = &t[0];

	t[3].edge[0].from  = &t[2];
	t[3].edge[1].to    = &t[1];
	t[3].edge[2].to    = &t[0];

	struct {int left; int right;} con[4][4];

	con [0][1].left  = 3;
	con [0][1].right = 2;

	con [0][2].left  = 1;
	con [0][2].right = 3;

	con [0][3].left  = 2;
	con [0][3].right = 1;

	con [1][0].left  = 2;
	con [1][0].right = 3;

	con [1][2].left  = 3;
	con [1][2].right = 0;

	con [1][3].left  = 0;
	con [1][3].right = 2;

	con [2][0].left  = 3;
	con [2][0].right = 1;

	con [2][1].left  = 0;
	con [2][1].right = 3;

	con [2][3].left  = 1;
	con [2][3].right = 0;

	con [3][0].left  = 1;
	con [3][0].right = 2;

	con [3][1].left  = 2;
	con [3][1].right = 0;

	con [3][2].left  = 0;
	con [3][2].right = 1;


	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			if (i==j) continue;
			for (k=0; k<2; k++)
			{
				if (t[j].edge[k].from != &t[i]) continue;
				if (t[j].edge[k].to != &t[j]) 
				{
					printf ("Error! bad edge setup!\n");
				}

				int tv = con[i][j].left;
				for (m=0;m<2; m++)
				{
					if (t[tv].edge[m].from == &t[j]) break;
				}
				t[j].edge[k].left = &t[tv].edge[m];

				tv = con[i][j].right;
				for (m=0;m<2; m++)
				{
					if (t[tv].edge[m].from == &t[j]) break;
				}
				t[j].edge[k].right = &t[tv].edge[m];
			}
		}
	}

	return t;
}

void talk_a_walk (Edge *e, int depth)
{
	if (0 >= depth) return;
	depth --;

	e->stats.cnt ++;

	talk_a_walk (e->left, depth);
	talk_a_walk (e->right, depth);
}

void show_stats (Vertex *t)
{
	int i, j;
	for (i=0; i<4; i++)
	{
		for (j=0; j<3; j++)
		{
			Edge *e = &t[i].edge[j];
			printf ("edge from %d to %d visited %d times\n", e->from->id, e->to->id, e->stats.cnt);
		}
	}
}

main ()
{
	Vertex *t = tetrahedron_new();
	Edge *e = &t->edge[0];

	talk_a_walk (e, 10);

	show_stats (t);
}
