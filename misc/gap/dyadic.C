
/*
 * dyadic.C
 *
 * Some exploratinos of polynomials generated from binary expansion
 *
 * Linas November 2004
 */

#include <stdio.h>

class Dyadic
{
	public:
		Dyadic (void);

		// assumes x = p/2^n 
		void SetFrac (unsigned int p, unsigned int order);

		// return alternating polynomial sum_n w^n (2b_n -1)
		double ToAlternatingPoly (double w);
	private:
		unsigned int ndigits;
		short bdigits[32];
};

Dyadic::Dyadic(void)
{
	ndigits = 0;
	int i;
	for (i=0; i<32; i++) {
		bdigits[i] = 0;
	}
}

 
void 
Dyadic :: SetFrac (unsigned int p, unsigned int order)
{
	ndigits = 0;
	if (0 == p) return;

	int sorder = (int) order;
	// make sure p isn't even, no powers of 2 in it
	while (p%2 == 0)
	{
		p /= 2;
		sorder --;
	}
	if (0 >= sorder) return;
	ndigits = sorder;

	unsigned int deno = 1<<(sorder-1);

	int i=0;
	while (1)
	{
		if (p <deno)
		{
			bdigits[i] = 0;
		}
		else
		{
			bdigits[i] = 1;
			p -= deno;
		}
		deno >>= 1;
		i++;
		if (0 == deno) break;
	}
}

double
Dyadic :: ToAlternatingPoly (double z)
{
	int i;
	double acc = 0.0;
	double zn = 1.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = 2*bdigits[i] - 1;
		double term = alt;
		term *= zn;
		acc += term;

		zn *= z;
	}

	return acc;
}


main (int argc, char * argv[])
{
	Dyadic dy;
	int order = 9;
	int nmax = 1<<order;

	double z = 0.3;
	int p;
	for (p=1; p<nmax; p+=2)
	{
		if (0 == p%2) continue;
		double x = ((double) p) / ((double) nmax);

		dy.SetFrac (p, order);
		double y = dy.ToAlternatingPoly (z);

		printf ("%d	%8.6g	%8.6g\n", p, x, y);
	}
}
