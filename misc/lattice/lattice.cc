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

	// multiply into this matrix
	void mult(const Uni&);

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
	a = 0; b = 1; c = -1; d = 0;
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

void Uni::mult(const Uni& o)
{
	long na = a * o.a + b * o.c;
	long nb = a * o.b + b * o.d;
	long nc = c * o.a + d * o.c;
	long nd = c * o.b + d * o.d;
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
	Uni u;

	double step = width / (double) sizex;
	// double yoff = height / (double) sizey;
	for (int k=0; k<itermax; k++)
	{
		unsigned long tk = 1<<k;
		for (unsigned long m=0; m<tk; m++)
		{
			u.make_uni(m, k);

			double x = re_center - 0.5*width + 0.5*step;
			// printf("sizex=%d sizey=%d center=(%g, %g) width=%g\n",
			//        sizex, sizey, re_center, im_center, width);
			for (int i = 0; i<sizex; i++)
			{
				double y = u.xform(x);
				y -= floor(y);

				x += step;
				if (0.0 > y) continue;
				if (1.0 <= y) continue;

				double ny = 0.5 + (y - im_center) / height;
				int j = floor(sizey * ny);
				j = (sizey-1) - j;
				if (j < 0)
				{
					printf("too small %g %g = %d %d\n", x-step, y, i, j);
					printf("uni= %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
					continue;
				}
				if (sizey <= j)
				{
					printf("too big %g %g = %d %d\n", x-step, y, i, j);
					continue;
				}
				glob[sizex*j+i] ++;

			}
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
