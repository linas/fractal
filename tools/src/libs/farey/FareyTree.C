
/*
 * FareyTree.C
 *
 * Walk the Farey Tree, returning the next Farey number in the tree.
 * This walk is a breadth-first walk, of course, alternating branches
 * until each depth has been explored (it does not go left-to-right
 * at each level).
 *
 * It uses memory in very intensive way.
 *
 * Created by Linas Vepstas October 2004
 */

#include <stdio.h>
#include <stdlib.h>
#include "FareyTree.h"

int
FareyIterator :: Recur (int np1, int dp1, int np2, int dp2, int *n, int *d, int depth)
{
	int nn, dd;
	nn = np1 + np2;
	dd = dp1 + dp2;

	/* Need to pack the state so that we don't waste memory too much */
	int pong = depth/4;
	int ping = depth%4;
	int pack = state[pong];
	pack >>= 2*ping;
	pack &= 0x3;
#define SET_STATE(ns) { int mask = 0x3 << (2*ping); state[pong] &= ~mask; state[pong] |= ns<<(2*ping); }

	if (0 == pack)
	{
		*n = nn;
		*d = dd;
		SET_STATE (1);
		return depth;
	}
	if (1 == pack)
	{
		SET_STATE (2);
		return Recur (np1, dp1, nn, dd, n, d, 3*depth+1);
	}
	if (2 == pack)
	{
		SET_STATE (1);
		return Recur (nn, dd, np2, dp2, n, d, 3*depth+2);
	}
	return depth;
}

void
FareyIterator :: GetNextFarey (int *n, int *d)
{
	int lvl = Recur (0,1, 1,1, n,d, 0);
	if (3*lvl+3>state_sz)
	{
		int newsz = (3*state_sz+7)/4;
		state = (char *) realloc (state, newsz*sizeof(char));
		if (!state)
		{
			fprintf (stderr, "going down, out of memory newsz=%d\n", newsz);
			exit (1);
		}
		int i;
		for (i=state_sz; i<newsz; i++)
		{
			state[i] = 0;
		}
		// printf ("duude realloc old=%d new=%d\n", state_sz, newsz);
		state_sz = 4*newsz;
	}
}

FareyIterator :: FareyIterator (void)
{
	#define FAREY_INITIAL_SIZE 11
	state = (char *) malloc (FAREY_INITIAL_SIZE *sizeof(char));

	int i;
	for (i=0; i<FAREY_INITIAL_SIZE; i++) 
	{
		state[i] = 0;
	}
	state_sz = 4*FAREY_INITIAL_SIZE;
}

FareyIterator :: ~FareyIterator ()
{
	free (state);
}

