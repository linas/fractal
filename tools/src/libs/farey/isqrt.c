/*
 * isqrt.c
 *
 * Return the integer square-root of an integer.
 *
 * Updates October 2016
 */

#include <math.h>

#include "isqrt.h"

/* ====================================================== */

long isqrt(long n)
{
	long res = 0;
	long bit = 1 << (8*sizeof(int)-2);
	while (bit > n) bit >>= 2;

	while (bit != 0)
	{
		long rpb = res + bit;
		if (n >= rpb)
		{
			n -= rpb;
			res = (res >> 1) + bit;
		}
		else  res >>= 1;
		bit >>= 2;
	}
	return res;
}

long integer_sqrt(long n)
{
	if (0 == n) return 0;
	long lo = 2 * integer_sqrt(n/4);
	long hi = lo+1;
	if (n < hi*hi) return lo;
	return hi;
}

static long integer_nth_helper(long n, long p, long tp)
{
	if (0 == n) return 0;
	long lo = 2 * integer_nth_helper(n/tp, p, tp);
	long hi = lo+1;

	long hipow = 1;
	for (long i=0; i<p; i++) hipow *= hi;
	if (n < hipow) return lo;
	return hi;
}

long integer_nth_root(long n, long p)
{
	return integer_nth_helper(n, p, 1<<p);
}

/* ====================================================== */

// #define TEST 1
#ifdef TEST

#include <stdio.h>

long test_isqrt (void)
{
	long have_error=0;
	long i;
	long nmax=10000;
	for (i=1; i<=nmax; i++)
	{
		long is = isqrt(i);
		if (is != integer_sqrt(i))
		{
			printf ("ERROR: in long sqrt function at n=%d\n", i);
			have_error ++;
		}
		if (is*is > i)
		{
			printf ("ERROR: in integer sqrt function at n=%d\n", i);
			have_error ++;
		}
	}

	if (0 == have_error)
	{
		printf ("PASS: tested isqrt function up to %d\n", nmax);
	}
	return have_error;
}

long test_isecond (void)
{
	long have_error=0;
	long i;
	long nmax=10000;
	for (i=1; i<=nmax; i++)
	{
		long is = isqrt(i);
		if (is != integer_nth_root(i, 2))
		{
			printf ("ERROR: in long 2nd root function at n=%d\n", i);
			have_error ++;
		}
		if (is*is > i)
		{
			printf ("ERROR: in integer 2nd root function at n=%d\n", i);
			have_error ++;
		}
	}

	if (0 == have_error)
	{
		printf ("PASS: tested 2nd root function up to %d\n", nmax);
	}
	return have_error;
}

long test_icube (void)
{
	long have_error=0;
	if (2 != integer_nth_root(8, 3)) have_error++;
	if (2 != integer_nth_root(26, 3)) have_error++;
	if (3 != integer_nth_root(27, 3)) have_error++;
	if (3 != integer_nth_root(63, 3)) have_error++;
	if (4 != integer_nth_root(64, 3)) have_error++;
	if (4 != integer_nth_root(124, 3)) have_error++;
	if (5 != integer_nth_root(125, 3)) have_error++;

	long nmax=10000;
	for (long i=1; i<=nmax; i++)
	{
		long is = integer_nth_root(i, 3);
		if (is*is*is > i)
		{
			printf ("ERROR: in integer cube root function at n=%d\n", i);
			have_error ++;
		}
	}

	if (0 == have_error)
	{
		printf ("PASS: tested cube root function up to %d\n", nmax);
	}
	return have_error;
}

int main()
{
	test_isqrt();
	test_isecond();
	test_icube();

	return 1;
}
#endif /* TEST */


/* --------------------------- END OF FILE ------------------------- */
