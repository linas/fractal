
/* crookedness sequence
 * WAYNE LEWIS AND PIOTR MINC
 * DRAWING THE PSEUDO-ARC
 *
 * Created Linas Vepstas November 2017
 */

#include <stdio.h>

int crook(int n)
{
	if (1 == n) return 1;
	if (2 == n) return 2;
	return 2*crook(n-1) + crook(n-2);
}

int main()
{
	for (int i=1; i<20; i++)
	{
		int cr = crook(i);
		double fcr = cr / ((double) (1<<(i-1)));
		printf("its %d %d %g\n", i, cr, fcr);
	}
}
