
/*
 * dyadic.C
 *
 * Some exploratinos of polynomials generated from binary expansion
 *
 * Linas November 2004
 * Linas December 2005
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

		// return polynomial sum_n b_n w^n 
		double ToCantorPoly (double w);

		// return alternating polynomial sum_n w^n (2b_n -1)
		double ToAlternatingPoly (double w);

		// return sum_n b_n w^n /n!
		double ToCantorExp (double w);

		// return sum_n (b_n XOR b_{n+1}) w^n
		double ToXOR (double w);

		// return sum_n (2b_n -1) n^-s
		double ToZetaPoly (double s);
		double ToZetaPolyC (double re_s, double im_s);

		// return product pi_n (b_n+1) 
		double ToProd (void);

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

// sum in an exponential-like fashion
double
Dyadic :: ToCantorExp (double z)
{
	int i;
	double acc = 0.0;
	double zn = 1.0;
	double fact = 1.0;

	for (i=0; i<ndigits; i++)
	{
		// short alt = bdigits[i];
		short alt = 2*bdigits[i] - 1;
		double term = alt;
		term *= zn;
		term /= fact;
		acc += term;

		zn *= z;
		fact *= i+1;
	}

	return acc;
}

double
Dyadic :: ToXOR (double z)
{
	int i;
	double acc = 0.0;
	double zn = 1.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = bdigits[i];
		if (alt == bdigits[i+1]) { alt=0; } else {alt=1;}
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
Dyadic :: ToProd (void)
{
	int i;
	double acc = 1.0;

	for (i=0; i<ndigits; i++)
	{
		short alt = bdigits[i] + 1;
		double term = alt;
		acc *= 0.5*term;
	}

	return acc;
}


double prod_acc (int nmax, int order)
{
	Dyadic dy;

	double acc =0.0;
	int k;
	for (k=1; k<nmax; k+=2)
	{
		int p = k;
		if (0 == p%2) continue;

		dy.SetFrac (p, order);
		double y = dy.ToProd();
		acc += y;
	}
	return acc;
}

int
main (int argc, char * argv[])
{
	Dyadic dy;
	int order = 10;
	int nmax = 1<<order;

	double zre = 1.0/3.0;
	double zim = 1.0/3.0;
	if (2 > argc) {
		printf ("Usage: %s <z value>\n", argv[0]); 
		exit(1);
	}
   zim = 0.0;
	zre = atof (argv[1]);
	printf ("# z=%g + i %g\n", zre, zim);

	double acc = 0.0;
	int k;
	for (k=1; k<nmax; k++)
	{
		int p = k;
		if (0 == p%2) continue;
		double x = ((double) p) / ((double) nmax);

#if 0
		int norder = order;
		while (0 == p%2)
		{
			p >>= 1;
			norder --;
		}
#endif

		dy.SetFrac (p, order);
		// double y = dy.ToAlternatingPoly (zre);
		// double y = dy.ToCantorPoly (zre);
		// double y = dy.ToCantorExp (zre);
		double y = dy.ToXOR (zre);
		// double y = dy.ToZetaPoly (zre);
		// double y = dy.ToZetaPolyC (zre, zim);

		acc += y;

#ifdef PASCAL_TAKAGI
		double y = dy.ToProd();
		double acc = prod_acc (p, order);
		double inv = prod_acc (nmax-p, order);
#endif

#ifdef COMB_STRUCTURE
		dy.SetFrac (p+1, order);
		y = dy.ToZetaPoly (z);
		dy.SetFrac (p-1, order);
		y -= dy.ToZetaPoly (z);
#endif

		printf ("%d	%8.6g	%8.6g	%8.6g\n", p, x, y, acc);
	}
	return 0;
}
