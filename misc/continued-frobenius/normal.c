
/*
 * normal.c
 *
 * Check to make sure that GKW is normal, or not.
 *
 * Linas Jan 2010
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ache.h"

// =================================================================

#define MS 100 // matrix dimension
#include "matrix.h"

// =================================================================

matrix  *gfp; // gkw 

void initialize (void)
{
	int i,j;

	printf ("sizeof matrix = %d\n", sizeof (matrix));
	
	gfp = (matrix *) malloc (sizeof (matrix));
	
	// initialize the frobenius-perron operator
	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			*gfp[i][j] = ache_mp(i,j);
		}
	}
}

// =================================================================

int
main ()
{
	initialize();
	return 0;
}
