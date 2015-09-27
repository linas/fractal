/*
 * lattice.cc
 *
 * Graph every line of the form (ax+b)/(cx+d) for SL(2,Z)
 * Do the same for the dyadic lines, as well.
 *
 * Linas Sept 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

struct Uni
{
	// m and k define a tree coordinate. m is a sequence of
	// left-right moves, encoded as bits; where 0==left move
	// 1==right move. k is the total number of moves, i.e. the
	// tree depth.
	unsigned long m;
	short k;

	// a,b,c,d define a (projective) unimodular matrix.
	long a,b,c,d;

	// Return S,T L R matrices
	void make_I(void);
	void make_S(void);
	void make_T(void);
	void make_L(void);
	void make_R(void);
	void make_P(void);
	void make_V(void);

	// multiply into this matrix
	void mult(const Uni&);

	// Compute the inverse
	void invert(void);

	// Given a tree coordinate, create the corresponding unimodular matrix
	void make_uni(unsigned long, short);

	// Mobius transform
	double xform(double);
};

void Uni::make_I(void)
{
	a = 1; b = 0; c = 0; d = 1;
}

void Uni::make_S(void)
{
	a = 0; b = -1; c = 1; d = 0;
}

void Uni::make_T(void)
{
	a = 1; b = 1; c = 0; d = 1;
}

void Uni::make_L(void)
{
	a = 1; b = 0; c = 1; d = 1;
}

// R is the same as T
void Uni::make_R(void)
{
	a = 1; b = 1; c = 0; d = 1;
}

void Uni::make_P(void)
{
	a = 0; b = -1; c = 1; d = 1;
}

// V is the same as S
void Uni::make_V(void)
{
	a = 0; b = -1; c = 1; d = 0;
}

void Uni::mult(const Uni& o)
{
	long na = a * o.a + b * o.c;
	long nb = a * o.b + b * o.d;
	long nc = c * o.a + d * o.c;
	long nd = c * o.b + d * o.d;
	a = na; b = nb; c = nc; d = nd;
}

void Uni::invert(void)
{
	long na = d;
	long nb = -b;
	long nc = -c;
	long nd = a;
	a = na; b = nb; c = nc; d = nd;
}

void Uni::make_uni(unsigned long em, short kay)
{
	Uni L, R;
	L.make_L();
	R.make_R();

	m = em; k = kay;
	make_I();

	// bit==0 == left move, else right move.
	for (short i=0; i<kay; i++)
	{
		if ((em & 0x1) == 0) mult(L);
		else mult(R);
		em >>= 1;
	}
}

double Uni::xform(double x)
{
	return (a*x+b) / (c*x+d);
}


void mobius_splat (
   float    *glob,
   int      sizex,
   int      sizey,
   double   re_center,
   double   im_center,
   double   width,
   double   height,
	Uni      u)
{
	double xstep = width / (double) sizex;
	double ystep = height / (double) sizey;

	// Run along x ...
	double x = re_center - 0.5*width + 0.5*xstep;
	// printf("sizex=%d sizey=%d center=(%g, %g) width=%g\n",
	//        sizex, sizey, re_center, im_center, width);
	for (int i = 0; i<sizex; i++)
	{
		double y = u.xform(x);
		y -= floor(y);

		x += xstep;
		if (0.0 > y) continue;
		if (1.0 <= y) continue;

		double ny = 0.5 + (y - im_center) / height;
		int j = floor(sizey * ny);
		j = (sizey-1) - j;
		if (j < 0)
		{
			printf("too small %g %g = %d %d\n", x-xstep, y, i, j);
			printf("uni= %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
			continue;
		}
		if (sizey <= j)
		{
			printf("too big %g %g = %d %d\n", x-xstep, y, i, j);
			continue;
		}
		glob[sizex*j+i] ++;
	}

	// Now fill in the other direction (anti-aliasing)
	u.invert();
	double y = re_center - 0.5*height + 0.5*ystep;
	// printf("sizey=%d sizex=%d center=(%g, %g) yidth=%g\n",
	//        sizey, sizex, re_center, im_center, yidth);
	for (int i = 0; i<sizey; i++)
	{
		double x = u.xform(y);
		x -= floor(x);

		y += ystep;
		if (0.0 > x) continue;
		if (1.0 <= x) continue;

		double nx = 0.5 + (x - im_center) / height;
		int j = floor(sizex * nx);
		j = (sizex-1) - j;
		if (j < 0)
		{
			printf("x too small %g %g = %d %d\n", y-ystep, x, i, j);
			printf("x uni= %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
			continue;
		}
		if (sizex <= j)
		{
			printf("x too big %g %g = %d %d\n", y-ystep, x, i, j);
			continue;
		}
		glob[sizey*i+j] ++;
	}
}

void MakeHisto (
   float    *glob,
   int      sizex,
   int      sizey,
   double   re_center,
   double   im_center,
   double   width,
   double   height,
   int      itermax,
   double   renorm)
{
	Uni I, V, P;
	I.make_I();
	V.make_V();
	P.make_P();

	for (int k=0; k<itermax; k++)
	{
		unsigned long tk = 1<<k;
		for (unsigned long m=0; m<tk; m++)
		{
			Uni a;
			a.make_uni(m, k);

			Uni u;
			u = I;
			u.mult(a);
			mobius_splat(glob, sizex, sizey, re_center, im_center, width, height, u);
			u.mult(V);
			mobius_splat(glob, sizex, sizey, re_center, im_center, width, height, u);

			u = P;
			u.mult(a);
			mobius_splat(glob, sizex, sizey, re_center, im_center, width, height, u);
			u.mult(V);
			mobius_splat(glob, sizex, sizey, re_center, im_center, width, height, u);
		}
	}
}


#if 0
int main(int argc, char* argv[])
{
	Uni u;

	for (int k=0; k<6; k++)
	{
		unsigned long tk = 1<<k;
		for (unsigned long m=0; m<tk; m++)
		{
			u.make_uni(m, k);

			// for (double x=0;
			printf("its %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
		}
	}
}
#endif
