/*
 * self.c
 *
 * Self-similarity explorer
 * Linas Vepstas Dec 2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void get_mob(int a1, int a2, int *a, int *b, int *c, int *d)
{

	*a = 1-a1*(a1-a2);
	*b = a1-a2;
	*c = -(a1-a2)*(1+a1*a2);
	*d = 1+a2*(a1-a2);
}

int det(int M[4][4])
{
	int d = 0;
	d += M[0][0] * M[1][1] * M[2][2] * M[3][3];
	d -= M[0][0] * M[1][1] * M[2][3] * M[3][2];
	d -= M[0][0] * M[1][2] * M[2][1] * M[3][3];
	d += M[0][0] * M[1][2] * M[2][3] * M[3][1];
	d += M[0][0] * M[1][3] * M[2][1] * M[3][2];
	d -= M[0][0] * M[1][3] * M[2][2] * M[3][1];

	d -= M[0][1] * M[1][0] * M[2][2] * M[3][3];
	d += M[0][1] * M[1][0] * M[2][3] * M[3][2];
	d += M[0][1] * M[1][2] * M[2][0] * M[3][3];
	d -= M[0][1] * M[1][2] * M[2][3] * M[3][0];
	d -= M[0][1] * M[1][3] * M[2][0] * M[3][2];
	d += M[0][1] * M[1][3] * M[2][2] * M[3][0];

	d += M[0][2] * M[1][1] * M[2][0] * M[3][3];
	d -= M[0][2] * M[1][1] * M[2][3] * M[3][0];
	d -= M[0][2] * M[1][0] * M[2][1] * M[3][3];
	d += M[0][2] * M[1][0] * M[2][3] * M[3][1];
	d += M[0][2] * M[1][3] * M[2][1] * M[3][0];
	d -= M[0][2] * M[1][3] * M[2][0] * M[3][1];

	d -= M[0][3] * M[1][0] * M[2][2] * M[3][1];
	d += M[0][3] * M[1][0] * M[2][1] * M[3][2];
	d += M[0][3] * M[1][2] * M[2][0] * M[3][1];
	d -= M[0][3] * M[1][2] * M[2][1] * M[3][0];
	d -= M[0][3] * M[1][1] * M[2][0] * M[3][2];
	d += M[0][3] * M[1][1] * M[2][2] * M[3][0];

	return d;
}

main (int argc, char * argv[])
{
	int arg = atoi(argv[1]);

	int a,b,c,d;
	int a1 = 1;
	int a2 = 3;
	a2=arg;
	get_mob (a1, a2, &a, &b, &c, &d);

	int aa,bb,cc,dd;

	for (a1=a2+1; a1<20; a1++)
	{
		get_mob (a1, a2, &aa, &bb, &cc, &dd);

		int m[4][4];
#define SIMILARITY
#ifdef SIMILARITY
		m[0][0] = a-aa;
		m[0][1] = -cc;
		m[0][2] = b;
		m[0][3] = 0;

		m[1][0] = -bb;
		m[1][1] = a-dd;
		m[1][2] = 0;
		m[1][3] = b;

		m[2][0] = c;
		m[2][1] = 0;
		m[2][2] = d-aa;
		m[2][3] = -cc;

		m[3][0] = 0;
		m[3][1] = c;
		m[3][2] = -bb;
		m[3][3] = d-dd;
#endif
#ifdef SIMPLE
		m[0][0] = a;
		m[0][1] = 0;
		m[0][2] = b;
		m[0][3] = 0;

		m[1][0] = 0;
		m[1][1] = a;
		m[1][2] = 0;
		m[1][3] = b;

		m[2][0] = c;
		m[2][1] = 0;
		m[2][2] = d;
		m[2][3] = 0;

		m[3][0] = 0;
		m[3][1] = c;
		m[3][2] = 0;
		m[3][3] = d;
#endif // SIMPLE

		int dt = det(m);

		int p4 = (a1-1)*(a1-1);
		p4 *= p4;
		printf("a1=%d det=%d p4=%d\n", a1, dt, p4);
	}
}
