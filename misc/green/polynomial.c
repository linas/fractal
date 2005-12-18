/*
 * Symbolic polynomial handling in c language
 * Very simple, only handles polynomials in one variable
 *
 * Linas Vepstas 16 December 2005
 */

#define MAXORDER 100

typedef double Poly[MAXORDER];

void 
poly_add (Poly *sum, Poly *a, Poly *b)
{
	int i;
	for (i=0; i<MAXORDER; i++)
	{
		(*sum)[i] = (*a)[i] + (*b)[i];
	}
}

