/*
 * Set up the matrix elements of the transfer operator 
 * (composition operator) of the Minkowski question mark function.
 * The basis are Haar wavelets.
 *
 * Linas Vepstas January 2006 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// return single-integer index corresponding to the
// discrete wavelet index (j,k)
#define HIDX(j,k)  ((1<<(j))+(k))

// This class is meant to return the farey number corresponding
// to the k/(2^j)'th dyadic fraction.
//
// Note that two competing coordinate systems are at play
// The "external" coord system are the values k/2^j 
// The "internal" coord system are the values (2k+1)/2^{j+1}
//
class SternBrocotTree
{
	public:
		SternBrocotTree (void);
		double GetFarey (int j, int k);

	private:
		void Enlarge (int j);
		void Fill (int jlo, int jhi);
		int max_j;
		int * numerators;
		int * denominators;
};

SternBrocotTree::SternBrocotTree (void)
{
	max_j = 4;
	int maxidx = 1<<max_j +3;
	numerators = (int *) malloc (maxidx * sizeof (int));
	denominators = (int *) malloc (maxidx * sizeof (int));
	Fill (0, max_j);
}

void SternBrocotTree::Enlarge (int j)
{
	j++;
	int maxidx = 1<<j +3;
	numerators = (int *) realloc (numerators, maxidx * sizeof (int));
	denominators = (int *) realloc (denominators, maxidx * sizeof (int));
	Fill (max_j, j);
	max_j = j;
}

void SternBrocotTree::Fill (int jlo, int jhi)
{
	int j,k;

	for (j=jlo; j<jhi;  j++)
	{
		numerators [HIDX(j,0)] = 1;
		denominators [HIDX(j, 0)] = j+2;

		if (0 == j) continue;

		numerators [HIDX(j,(1<<j)-1)] = j+1;
		denominators [HIDX(j, (1<<j)-1)] = j+2;
		for (k=1; k<(1<<j)-1; k++)
		{
			int kl, kr, jl, jr; 
			kl = k;
			kr = k+1;
			jl = j-1;
			jr = j-1;
			while (kl%2 == 0) {kl>>= 1; jl -= 1; } 
			while (kr%2 == 0) {kr>>= 1; jr -= 1; } 
			kl = (kl-1)>>1;
			kr = (kr-1)>>1;
			
			numerators [HIDX(j,k)] = numerators [HIDX(jl,kl)] + numerators [HIDX(jr,kr)];
			denominators [HIDX(j,k)] = denominators [HIDX(jl,kl)] + denominators [HIDX(jr,kr)];
		}
	}

#if DEBUG
	for (j=jlo; j<jhi;  j++)
	{
		for (k=0; k<(1<<j); k++)
		{
			int i = HIDX(j,k);
			printf ("really its idx=%d (%d/%d) == (%d/%d)\n", HIDX(j,k), 2*k+1, 1<<(j+1), 
			        numerators [HIDX(j,k)], denominators [HIDX(j,k)]);
		}
		printf ("\n");
	}
#endif
}

double SternBrocotTree::GetFarey (int j, int k)
{
	if (0 == k) return 0.0;
	j--;

	while (k%2 == 0) { k>>= 1; j -= 1; } 
	k = (k-1)/2;

	if (0 > j) return 1.0;

	if (j >= max_j) Enlarge (j);

	return ((double)  numerators [HIDX(j,k)]) / ((double) denominators [HIDX(j,k)]);
}

// ==========================================================================

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

inline double
haar_domain_midpoint (int j, int k)
{
	return ((double) 2*k+1) / ((double) (1<<(j+1)));
}

// return overlap of two square bumps
static inline double 
overlap (double alo, double ahi, double blo, double bhi)
{
	double bot=alo;
	double top = ahi;
	if (bot < blo) bot = blo;
	if (bhi < top) top = bhi;
	if (top <= bot) return 0.0;
	return top-bot;
} 

// return the matrix element of the question mark 
// composition operator (== transfer operator).
double 
q_oper_elt (int j, int k, int l, int m)
{
	static SternBrocotTree t;

	// There are two non-intersection cases
	double haar_hi = haar_domain_max (l,m);
	double farey_lo = t.GetFarey (j, k);
	if (haar_hi <= farey_lo) return 0.0;

	double haar_lo = haar_domain_min (l,m);
	double farey_hi = t.GetFarey (j, k+1);
	if (haar_lo >= farey_hi) return 0.0;

	// There are four total-overlap cases
	double haar_mid = haar_domain_midpoint (l,m);
	double farey_mid = t.GetFarey (j+1, 2*k+1);

	if (haar_lo  <= farey_lo && farey_hi <= haar_mid) return 0.0;
	if (haar_mid <= farey_lo && farey_hi <= haar_hi) return 0.0;
	if (farey_lo  <= haar_lo && haar_hi <= farey_mid) return 0.0;
	if (farey_mid <= haar_lo && haar_hi <= farey_hi) return 0.0;

	// printf ("haar= %g %g %g\n", haar_lo, haar_mid, haar_hi);
	// printf ("farey= %g %g %g\n", farey_lo, farey_mid, farey_hi);

	double elt = 0.0;
	elt += overlap (haar_lo, haar_mid, farey_lo,  farey_mid);
	elt -= overlap (haar_lo, haar_mid, farey_mid, farey_hi);
	elt -= overlap (haar_mid, haar_hi, farey_lo,  farey_mid);
	elt += overlap (haar_mid, haar_hi, farey_mid, farey_hi);
	int p = j+l;
	if (p%2 == 1) elt *= M_SQRT2;
	elt *= 1 << (p/2);
	return elt;
}

int test (void)
{
	double x = q_oper_elt (1,1,2,0);
	double y = q_oper_elt (1,1,2,1);
	double z = q_oper_elt (1,1,2,2);
	double w = q_oper_elt (1,1,2,3);

	// should be 0.0 except for third one which is sqrt(8)/6
	printf ("its %g %g %g %g \n", x,y, z, w);

	int j=4;
	int k;
	for (k=0; k<(1<<j); k++)
	{
		double x = q_oper_elt (1,0,j,k);
		printf ("its %d %g\n", k, x);
	}
}

// ==========================================================================

void make_oper (int dim)
{
	dim ++;

	double *mat = (double *) malloc ((dim+1)*(dim+1)*sizeof (double));

	int j,k,l,m;

	for (j=0; j<dim; j++)
	{
		for (k=0; k<(1<<j); k++)
		{
			printf ("jk=%d %d == ",j,k); 
			for (l=0; l<dim; l++)
			{
				for (m=0; m<(1<<l); m++)
				{
					double x = q_oper_elt (j,k,l,m);
					printf ("%5.2g ", x);
				}
			}
			printf ("\n");
		}
	}
}

main ()
{
	make_oper (4);
return 1;
	int j=8;
	int k;
	for (k=0; k<(1<<j); k++)
	{
		double x = q_oper_elt (j,k, 3,0);
		printf ("its %d %g\n", k, x);
	}
}
