/*
 * Semi-Symbolic polynomial handling in c language
 * Very simple, only handles polynomials in one variable
 *
 * Linas Vepstas 16 December 2005
 */

#define MAXORDER 100

typedef double Poly[MAXORDER];

void 
poly_clear (Poly *b)
{
	int i;
	for (i=0; i<MAXORDER; i++)
	{
		(*b)[i] = 0.0;
	}
}

void 
poly_add (Poly *sum, Poly *a, Poly *b)
{
	int i;
	for (i=0; i<MAXORDER; i++)
	{
		(*sum)[i] = (*a)[i] + (*b)[i];
	}
}

void 
poly_mult (Poly *prod, Poly* a, Poly *b)
{
	int i,j;
	
	for (i=0; i<MAXORDER; i++)
	{
		double acc = 0.0;
		for (j=0; j<=i; j++)
		{
			acc += (*a)[j] * (*b)[i-j];
		}
		(*prod)[i] = acc;
	}
}

void 
poly_differntiate (Poly *deriv, Poly *in)
{
	int i;
	for (i=0; i<MAXORDER-1; i++)
	{
		(*deriv)[i] = (i+1) * (*in)[i+1];
	}
}

double 
poly_eval (Poly *a, double val)
{
	int i;
	double acc = 0.0;
	double term = 1.0;
	for (i=0; i<MAXORDER; i++)
	{
		acc += term * (*a)[i];
		term *= val;
	}

	return acc;
}

#define MAXTERMS 100
static Poly gaussian_derivatives[MAXTERMS];

void 
setup_gaussians (double sigma, double mu)
{
	int i;
	for(i=0; i<MAXTERMS; i++)
	{
		poly_clear (&gaussian_derivatives[i]);
	}
	gaussian_derivatives[0][0] = 1.0;
	gaussian_derivatives[1][0] = mu/(sigma*sigma);
	gaussian_derivatives[1][1] = -mu/(sigma*sigma);
	

}

