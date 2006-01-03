/*
 * Set up the matrix elements of the transfer operator 
 * (composition operator) of the Minkowski question mark function.
 * The basis are Haar wavelets.
 *
 * Linas Vepstas January 2006 
 */

// return single-integer index corresponding to the
// discrete wavelet index (j,k)
#DEFINE HIDX(j,k)  ((1<<((j)-1))+(k))

// Return the minimum value below which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_min (int j, int k)
{
	return ((double) k) / ((double) (1<<(j))
}

// Return the maximum value above which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_max (int j, int k)
{
	return ((double) k+1) / ((double) (1<<(j))
}

class SternBrocotTree
{
	public:
		SternBrocotTree (void);
		double GetFarey (int j, int k);

	private:
		void Fill (int jlo, int jhi);
		int maxidx;
		int max_j;
		int * numerators;
		int * denominators;
};

SternBrocotTree::SternBrocotTree (void);
{
	max_j = 10;
	maxidx = 1<<max_j;
	numerators = (int *) malloc (maxidx * sizeof (int));
	denominators = (int *) malloc (maxidx * sizeof (int));
	Fill (0, max_j);
}

void SternBrocotTree::Fill (int jlo, int jhi)
{
	int j,k;

	for (j=jlo; j<=jhi;  j++)
	{
		for (k=0; k<(1<<j); k++)
		{
		}
	}
}

inline double
farey_domain_midpoint (int j, int k)
{
}

main ()
{
}
