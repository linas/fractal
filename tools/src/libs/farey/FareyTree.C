
/*
 * FareyTree.h
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
class FareyIterator
{
	public:
		FareyIterator (void);
		~FareyIterator ();
		void GetNextFarey (int *num, int *denom);
	protected:
		int Recur (int np1, int dp1, int np2, int dp2, int *n, int *d, int depth);
	private:
		char *state;
		int state_sz;
};

int
FareyIterator :: Recur (int np1, int dp1, int np2, int dp2, int *n, int *d, int depth)
{
	int nn, dd;
	nn = np1 + np2;
	dd = dp1 + dp2;
	if (0 == state[depth])
	{
		*n = nn;
		*d = dd;
		state[depth] = 1;
		return depth;
	}
	if (1 == state[depth])
	{
		
		state[depth] = 2;
		return Recur (np1, dp1, nn, dd, n, d, 3*depth+1);
	}
	if (2 == state[depth])
	{
		state[depth] = 1;
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
		int newsz = 3*state_sz+3;
		state = (int *) realloc (state, newsz*sizeof(int));
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
		state_sz = newsz;
	}
}

FareyIterator :: FareyIterator (void)
{
	#define FAREY_INITIAL_SIZE 3
	state = (int *) malloc (FAREY_INITIAL_SIZE *sizeof(int));

	int i;
	for (i=0; i<FAREY_INITIAL_SIZE; i++) 
	{
		state[i] = 0;
	}
	state_sz = FAREY_INITIAL_SIZE;
}

FareyIterator :: ~FareyIterator ()
{
	free (state);
}

