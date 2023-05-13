
/* density.C:
 *
 * Working over the set of rationals, determine the lengths of the corresponding
 * continued fractions.  I.E. What fraction of rationals have length 1,2,3 etc?
 * What is the mean length? (Note: the length of the continued fraction for a 
 * rational is always finitie).
 *
 * Linas Vepstas February 2004
 */


#include "Farey.h"

/* Implement Euclid's algorithm for greatest common divisor */
/* http://www.jimloy.com/number/euclids.htm */
int gcd (int dividend, int divisor)
{

	int remainder;

	while (1)
	{
		remainder = dividend % divisor;
		if (0==remainder) return divisor;
		dividend = divisor;
		divisor = remainder;
	}
	return 1;
}   

/* the lcm of a and b is equal to ab divided by the gcd of a and b */

void count_regular (ContinuedFraction *cf, int deno, int *distro)
{
	int num;
	for (num=1; num<deno; num++)
	{
		/* work only with relatively prime things */
		if (1 != gcd (deno,num)) continue;
		cf->SetRatio (num,deno);

		int nt = cf->GetNumTerms();
		distro[nt] ++;	
	}
}

#define NARR 200
void print_distro (int *distro)
{
	int cnt = 0;
	int i;
	for (i=0; i<NARR; i++)
	{
		cnt += distro[i];
	}
	for (i=2; i<16; i+=2)
	{
		printf ("%8.6f\t", ((double) distro[i])/((double) cnt));
	}
	printf ("\n");
}

void count (void)
{
	int distribution[NARR];

	int i;
	for (i=0; i<NARR; i++)
	{
		distribution[i] = 0;
	}
	
	ContinuedFraction *cf = new ContinuedFraction;
#if 0
	cf->SetRatio (2,3);
	printf ("its %d\n", cf->GetNumTerms());
	cf->Print ();
	return;
#endif
	
	int deno;
	for (deno=2; deno<1000000; deno++)
	{
		int cnt=0;
		int nterms=0;
		count_regular (cf, deno, distribution);

		if (0==deno%323) 
		{
			printf ("deno=%d ", deno);
			print_distro (distribution);
		}
	}
}

main ()
{
	count();
}
