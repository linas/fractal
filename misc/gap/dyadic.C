
/*
 * dyadic.C
 *
 * Some exploratinos of polynomials generated from binary expansion
 *
 * Linas November 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

class Dyadic
{
	public:
		Dyadic (void);

		// assumes x = p/2^n where n==order 
		// if p isn't odd, then this removes the common factors,
		// decrmenting the order as needed.
		// this function computes the binary expansion of p 
		// that is, gets the binary digits of p
		void SetFrac (unsigned int p, unsigned int order);

		// return alternating polynomial sum_n w^n (2b_n -1)
		double ToAlternatingPoly (double w);
		double ToCantorPoly (double w);

		// return sum_n (2b_n -1) n^-s
		double ToZetaPoly (double s);
		double ToZetaPolyC (double re_s, double im_s);

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


double
Dyadic :: ToZetaPoly (double s)
{
	int i;
	double acc = 0.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = 2*bdigits[i] - 1;
		double term = alt;

		double en = i+1;
		double z = pow (en, -s);
		term *= z;
		acc += term;
	}

	return acc;
}

double
Dyadic :: ToZetaPolyC (double s_re, double s_im)
{
	int i;
	double re_acc = 0.0;
	double im_acc = 0.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = 2*bdigits[i] - 1;
		double term = alt;

		double en = i+1;
		double z = pow (en, -s_re);
		term *= z;

		double arg = s_im * log (en);
		double co = cos (arg);
		double si = -sin (arg);
		re_acc += term*co;
		im_acc += term*si;
	}

	return re_acc;
}


double
Dyadic :: ToCantorPoly (double z)
{
	int i;
	double acc = 0.0;
	double zn = 1.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = bdigits[i];
		double term = alt;
		term *= zn;
		acc += term;

		zn *= z;
	}

	return acc;
}


int
main (int argc, char * argv[])
{
	Dyadic dy;
	int order = 12;
	int nmax = 1<<order;

	double zre = 1.0/3.0;
	double zim = 1.0/3.0;
	if (2 > argc) {
		printf ("Usage: %s <z value>\n", argv[0]); 
		exit(1);
	}
   zre = 1.3;
	zim = atof (argv[1]);
	printf ("# z=%g + i %g\n", zre, zim);

	int p;
	for (p=1; p<nmax; p+=2)
	{
		if (0 == p%2) continue;
		double x = ((double) p) / ((double) nmax);

		dy.SetFrac (p, order);
		// double y = dy.ToAlternatingPoly (z);
		// double y = dy.ToCantorPoly (z);
		double y = dy.ToZetaPolyC (zre, zim);

#ifdef COMB_STRUCTURE
		dy.SetFrac (p+1, order);
		y = dy.ToZetaPoly (z);
		dy.SetFrac (p-1, order);
		y -= dy.ToZetaPoly (z);
#endif

		printf ("%d	%8.6g	%8.6g\n", p, x, y);
	}
	return 0;
}
