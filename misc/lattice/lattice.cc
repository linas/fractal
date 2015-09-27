
/*
 * lattice.cc
 *
 * Graph every line of the form (ax+b)/(cx+d) for SL(2,Z)
 * Do the same for the dyadic lines, as well.
 *
 * Linas Sept 2015
 */

#include <stdio.h>

struct Uni
{
	// m and k define a dyadic rational m/2^k
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
	void make_uni(unsigned long, short);
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

int main(int argc, char* argv[])
{
	Uni u;

	u.make_uni(0, 2);
	printf("its %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
	u.make_uni(1, 2);
	printf("its %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
	u.make_uni(2, 2);
	printf("its %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
	u.make_uni(3, 2);
	printf("its %ld %ld %ld %ld\n", u.a, u.b, u.c, u.d);
}
