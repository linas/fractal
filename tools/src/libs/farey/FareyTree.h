
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

