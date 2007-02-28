
/*
 * treemap.c
 *
 * Map a binary tree onto a polyhedron
 *
 * More or less equivalent to a Markov chain.
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

		for (j=0; j<3; j++)
		{
			t[i].edge[j].stats.cnt = 0;
			t[i].edge[j].stats.sum = 0.0;
			t[i].edge[j].stats.sumsq = 0.0;
			t[i].edge[0].from = &t[i];
			t[i].edge[1].from = &t[i];
			t[i].edge[2].from = &t[i];
		}

		k = 0;
		for (j=0; j<4; j++)
		{
			if (j == i) continue;
			t[i].edge[k].to = &t[j];
			k++;
		}
	}

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
			for (k=0; k<3; k++)
			{
	// printf ("from %d to %d maybe by =%d->%d\n", i, j, t[i].edge[k].from->id, t[i].edge[k].to->id); 
				if (t[i].edge[k].to == &t[j]) break;
			}

			int tv = con[i][j].left;
			for (m=0;m<3; m++)
			{
				if (t[j].edge[m].to == &t[tv]) break;
			}
			t[i].edge[k].left = &t[j].edge[m];

// printf ("from %d to %d conected by =%d (%d->%d) want left  %d (%d->%d)\n", i,j, k, 
// t[i].edge[k].from->id, t[i].edge[k].to->id, tv,
// t[j].edge[m].from->id, t[j].edge[m].to->id);

			tv = con[i][j].right;
			for (m=0;m<3; m++)
			{
				if (t[j].edge[m].to == &t[tv]) break;
			}
			t[i].edge[k].right = &t[j].edge[m];

// printf ("from %d to %d conected by =%d (%d->%d) want right %d (%d->%d)\n\n", i,j, k, 
// t[i].edge[k].from->id, t[i].edge[k].to->id, tv,
// t[j].edge[m].from->id, t[j].edge[m].to->id);

		}
	}

	return t;
}

void talk_a_midpoint_walk (Edge *e, int depth, double left, double right)
{
	if (0 >= depth) return;
	depth --;

	double mid = 0.5*(left+right);

	e->stats.cnt ++;

	Vertex * v = e->to;
	v->stats.cnt ++;
	v->stats.sum += mid;
	v->stats.sumsq += mid*mid;

	// printf ("depth=%d came to %d\n", depth, e->to->id);

	talk_a_midpoint_walk (e->left, depth, left, mid);
	talk_a_midpoint_walk (e->right, depth, mid, right);
}

void talk_a_walk (Edge *e, int depth, double prob, double val)
{
	if (0 >= depth) return;
	depth --;

	e->stats.cnt ++;

	Vertex * v = e->to;
	v->stats.cnt ++;
	v->stats.sum += val;
	v->stats.sumsq += val*val;

	talk_a_walk (e->left, depth, prob, prob*val);
	talk_a_walk (e->right, depth, prob, (1.0-prob)*val);
}

void show_stats (Vertex *t, int depth)
{
	int cnt = 1<<depth;
	printf ("depth=%d cnt=%d\n", depth, cnt);

	double frac = 1.0 / ((double) cnt);
	double flen = 1.0/ ((double) depth);
	
	int i, j;
	int prev = 0;
	for (i=0; i<4; i++)
	{
		Vertex *v = &t[i];
		printf ("vertex %d visited %d times ", v->id, v->stats.cnt-prev);
		if (i==0) prev = v->stats.cnt;

		printf ("sum=%g ", flen * v->stats.sum);
		printf ("avg=%g ", frac * v->stats.sum);
		printf ("\n");
	}
	prev = 0;
	for (i=0; i<4; i++)
	{
		for (j=0; j<3; j++)
		{
			Edge *e = &t[i].edge[j];
			// printf ("edge from %d to %d visited %d times\n", e->from->id, e->to->id, e->stats.cnt-prev);
			if ((i==0) && (j==0)) prev = e->stats.cnt;
		}
	}
	printf ("\n");
}

int
main (int argc, char * argv[])
{
	if (argc<2)
	{
		fprintf (stderr, "Usage: %s <depth>\n", argv[0]);
		exit (1);
	}

	int depth = atoi (argv[1]);

	Vertex *t = tetrahedron_new();
	Edge *e = &t[3].edge[0];

	double left = 0.1;
	double right = 1.0;

	talk_a_walk (e, depth, left, right);
	show_stats (t, depth);
}
