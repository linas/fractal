/*
 * Set up the matrix elements of the transfer operator 
 * (composition operator) of the Minkowski question mark function.
 * The basis are Haar wavelets.
 *
 * Linas Vepstas January 2006 
 */

#include <stdio.h>
#include <stdlib.h>

// return single-integer index corresponding to the
// discrete wavelet index (j,k)
#define HIDX(j,k)  ((1<<(j))+(k))

// Return the minimum value below which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_min (int j, int k)
{
	return ((double) k) / ((double) (1<<(j)));
}

// Return the maximum value above which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_max (int j, int k)
{
	return ((double) k+1) / ((double) (1<<(j)));
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

SternBrocotTree::SternBrocotTree (void)
{
	max_j = 4;
	maxidx = 1<<(max_j+1);
	numerators = (int *) malloc (maxidx * sizeof (int));
	denominators = (int *) malloc (maxidx * sizeof (int));
	Fill (0, max_j);
}

void SternBrocotTree::Fill (int jlo, int jhi)
{
	int j,k;

	for (j=jlo; j<=jhi;  j++)
	{
		numerators [HIDX(j,0)] = 1;
		denominators [HIDX(j, 0)] = j+2;

printf ("wtf idx=%d (%d/%d)\n", HIDX(j,0), numerators [HIDX(j,0)], denominators [HIDX(j,0)]);
		if (0 == j) continue;

		numerators [HIDX(j,(1<<j)-1)] = j+1;
		denominators [HIDX(j, (1<<j)-1)] = j+2;
k=(1<<j)-1;
printf ("other idx=%d (%d/%d)\n", HIDX(j,k), numerators [HIDX(j,k)], denominators [HIDX(j,k)]);
		for (k=1; k<(1<<j)-1; k++)
		{
			int kl, kr, jl, jr; 
			kl = k;
			kr = k+1;
			jl = j-1;
			jr = j-1;
printf ("duude start left idx=%d (jl, kl)=(%d %d) lidx=%d \n", HIDX(j,k), jl, kl,  HIDX(jl,kl));
			while (kl%2 == 0) {kl>>= 1; jl -= 1; } 
			while (kr%2 == 0) {kr>>= 1; jr -= 1; } 
			kl = (kl-1)>>1;
			kr = (kr-1)>>1;
			
printf ("duude left idx=%d (jl, kl)=(%d %d) lidx=%d (%d/%d) \n", HIDX(j,k), jl, kl,  HIDX(jl,kl), numerators [HIDX(jl,kl)], denominators [HIDX(jl,kl)] );
printf ("duude rigt idx=%d (jr, kr)=(%d %d) ridx=%d (%d/%d) \n", HIDX(j,k), jr, kr,  HIDX(jr,kr), numerators [HIDX(jr,kr)], denominators [HIDX(jr,kr)] );
			numerators [HIDX(j,k)] = numerators [HIDX(jl,kl)] + numerators [HIDX(jr,kr)];
			denominators [HIDX(j,k)] = denominators [HIDX(jl,kl)] + denominators [HIDX(jr,kr)];

printf ("dude its idx=%d (%d/%d) == (%d/%d)\n\n", HIDX(j,k), 2*k+1, 1<<(j+1), numerators [HIDX(j,k)], denominators [HIDX(j,k)]);
		}
	}

printf ("\n");
	for (j=jlo; j<=jhi;  j++)
	{
		for (k=0; k<(1<<j); k++)
		{
			int i = HIDX(j,k);
printf ("really its idx=%d (%d/%d) == (%d/%d)\n", HIDX(j,k), 2*k+1, 1<<(j+1), numerators [HIDX(j,k)], denominators [HIDX(j,k)]);
		}
printf ("\n");
	}
	
}

inline double
farey_domain_midpoint (int j, int k)
{
}

main ()
{
	SternBrocotTree t;
}
