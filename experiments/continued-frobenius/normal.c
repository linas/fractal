
/*
 * normal.c
 *
 * try to see how non-normal the GKW is .. 
 * (well, its very non-normal....
 *
 * Linas Jan 2010
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ache.h"

// =================================================================

#define MS 10 // matrix dimension
#include "matrix.h"

// =================================================================

matrix  *gfp; // gkw 
matrix  *transp; // gkw 
matrix  *pa, *pb, *norm; // gkw 

void initialize (void)
{
	int i,j;
	
	printf ("sizeof matrix = %d\n", sizeof (matrix));

	gfp = (matrix *) malloc (sizeof (matrix));
	transp = (matrix *) malloc (sizeof (matrix));
	pa = (matrix *) malloc (sizeof (matrix));
	pb = (matrix *) malloc (sizeof (matrix));
	norm = (matrix *) malloc (sizeof (matrix));
	
	// initialize the frobenius-perron operator
	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
			(*gfp)[i][j] = ache_mp(i,j);
			(*transp)[i][j] = (*gfp)[j][i];
if (i<5 && j<5) printf("M[%d][%d] = %Lf\n", i, j, (*gfp)[i][j]);
		}
	}
	printf ("done init\n");
}

// =================================================================

int
main ()
{
	int i,j;

	initialize();
	multmatrix (pa, gfp, transp);
	multmatrix (pb, transp, gfp);
	scalematrix(pb, -1.0);
	addmatrix (norm, pa, pb);

	for (i=0; i<MS; i++)
	{
		for (j=0; j<MS; j++)
		{
if (i<5 && j<5) printf("norm[%d][%d] = %Lf\n", i, j, (*norm)[i][j]);
		}
	}
	return 0;
}
