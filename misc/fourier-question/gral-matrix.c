/*
 * Compute matrix elements of fourier of minkowski question mark 
 * by means of brute-force integration. This seems to be a viable
 * quick-n-dirty way of getting these.
 *
 * Linas June 2008
 */

#define M_PIl    3.1415926535897932384626433832795029L  /* pi */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int arrsz;
static long double *sin_arr;
static long double *cos_arr;
static long double delta;

void init_sine_cosine_arrays(int size)
{
	int i;
	arrsz = size;
	sin_arr = (long double *) malloc(arrsz * sizeof(long double));
	cos_arr = (long double *) malloc(arrsz * sizeof(long double));

	delta = 2.0L*M_PIl;
	delta /= (long double) arrsz;
	long double theta = 0.0L;
	for (i=0; i<arrsz; i++)
	{
		theta = ((long double)i)*delta;
		sin_arr[i] = sinl(theta);
		cos_arr[i] = cosl(theta);
	}
	printf ("# done with init\n");
}

void make_elt(int m, int n, long double *pre, long double *pim)
{
	long double start, end, th, off, g;
	long double re = 0.0L;
	long double im = 0.0L;
	int i, j;

	off = 16.0L * M_PIl*M_PIl * m * n;
	start = 2.0L*M_PIl * (m+n);
	end = M_PIl * (n+4*m);
	i = 0;
	th = start;
	while (th < end)
	{
		th = start + ((long double)i)*delta;
		g = th*th - off;
		g = th / sqrtl(g);

		j = i % arrsz;
		re += g* cos_arr[j];
		im -= g* sin_arr[j];

		i++;
	}

	re *= delta;
	im *= delta;

	re *= -1.0L / (4.0L * M_PIl * m);
	im *= -1.0L / (4.0L * M_PIl * m);

	if (n%2)
	{
		im += -1.0L / (2.0L * M_PIl * m);
	}

	*pre = re;
	*pim = im;
}

main()
{
	init_sine_cosine_arrays(4123123);

	long double re = 0.0L;
	long double im = 0.0L;
	make_elt(2,5, &re, &im);

	printf("its %Lg +i Lg\n", re, im);
}
