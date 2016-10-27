/*
 * moebius.c
 *
 * Return the Moebius function of an integer.
 * Not intended for large integers, works only for small integers,
 * due to poor-man's factorization algo.
 *
 * Linas Vepstas January 2005
 * Updates July 2006
 * Updates November 2014
 * Updates October 2016
 */

#include <math.h>
#include <malloc.h>

#include "gcf.h"
#include "moebius.h"
#include "cache.h"

/* ====================================================== */

int isqrt(int n)
{
	int res = 0;
	int bit = 1 << (8*sizeof(int)-2);
	while (bit > n) bit >>= 2;

	while (bit != 0)
	{
		int rpb = res + bit;
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

int integer_sqrt(int n)
{
	if (0 == n) return 0;
	int lo = 2 * integer_sqrt(n/4);
	int hi = lo+1;
	if (n < hi*hi) return lo;
	return hi;
}

static int integer_nth_helper(int n, int p, int tp)
{
	if (0 == n) return 0;
	int lo = 2 * integer_nth_helper(n/tp, p, tp);
	int hi = lo+1;

	int hipow = 1;
	for (int i=0; i<p; i++) hipow *= hi;
	if (n < hipow) return lo;
	return hi;
}

int integer_nth_root(int n, int p)
{
	return integer_nth_helper(n, p, 1<<p);
}

/* ====================================================== */

static unsigned int *sieve = NULL;
static unsigned int sieve_size = 0;
static unsigned int sieve_max = 0;

#define INIT_PRIME_SIEVE(N) \
	if (!sieve || sieve[sieve_max]*sieve[sieve_max] <(N)) {\
		init_prime_sieve(N); \
	}

/* Initialize and fill in a prime-number sieve.
 * Handles primes up to 4 billion (2^32)
 * long long int should be a 64-bit number
 *
 * XXX Ths should be converted to use the code in prime.h and prime.c,
 * but for now, its here, soe as to maximize copiler optimization.
 */
static void
init_prime_sieve (long long int prod)
{
	unsigned int n, j;
	unsigned int nstart;
	unsigned int pos;

	if (sieve)
	{
		long long int ss = sieve[sieve_max];
		if (ss*ss > prod) return;
	}

	unsigned int max = 1000.0+sqrt (prod);

	if (!sieve)
	{
		sieve = (unsigned int *) malloc (8192*sizeof (unsigned int));
		sieve_size = 8192;
		sieve_max = 2;
		sieve[0] = 2;
		sieve[1] = 3;
		sieve[2] = 5;
	}
	pos = sieve_max+1;
	nstart = sieve[sieve_max] + 2;

	/* dumb algo, brute-force test all odd numbers against known primes */
	for (n=nstart; n<=max; n+=2)
	{
		for (j=1; ; j++)
		{
			int p = sieve[j];
			if (0 == n%p)
			{
				break;
			}
			if (p*p > n)
			{
				/* If we got to here, n must be prime; save it, move on. */
				sieve[pos] = n;
				pos ++;
				if (pos >= sieve_size)
				{
					sieve_size += 8192;
					sieve = (unsigned int *)realloc (sieve, sieve_size * sizeof (unsigned int));
				}
				break;
			}
		}
	}
	sieve_max = pos-1;

#if 0
	for (j=0; j<pos; j++)
	{
		printf ("its %d %d\n", j, sieve[j]);
	}
#endif
}

/* ====================================================== */
/** Compute the divisor arithmetic function
 *  Returns the number of divisors of n.
 *  Almost a tail-recursive algorithm.
 */

static int divisor_helper (long long int n, int last_prime_checked)
{
	int ip;

	if (1==n) return 1;
	if (2==n) return 2;

	INIT_PRIME_SIEVE(n);

	for (ip=last_prime_checked+1; ; ip++)
	{
		long long int d = sieve[ip];

		if (d*d >n) break;
		if (n%d) continue;

		int acc = 2;
		n /=d;
		while (0 == n%d)
		{
			n /=d;
			acc ++;
		}
		if (1==n) return acc;
		return acc * divisor_helper (n, ip);
	}

	return 2;
}

int divisor (long long int n)
{
	return divisor_helper (n, -1);
}

/** 
 * Sigma arithmetic series, equals divisor arith series for a=0
 * Computes the divisors of n, raises each to the a'th power, and
 * returns thier sum.
 *
 * A slow, simple-minded algo.
 */
int sigma (int n, int a)
{
	int acc = 0;
	int d;

	int ns = isqrt(n);
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;

		int dp = 1;
		for (int ia=0; ia<a; ia++) dp *= d;
		acc += dp;
	}

	int dp = 1;
	for (int ia=0; ia<a; ia++) dp *= n;
	acc += dp;

	return acc;
}

/* same as above, but for float-point power */
long double sigmaf (int n, long double a)
{
	long double acc = 0;
	int d;

	int ns = isqrt(n);
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;
		acc += powl(d, a);
	}
	acc += powl(n, a);

	return acc;
}

/* same as above, but includes a log */
long double sigmalog (int n, long double a)
{
	long double acc = 0;
	int d;

	int ns = isqrt(n);
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;
		acc += powl(d, a) * logl(d);
	}
	acc += powl(n, a) * logl(n);

	return acc;
}

/* ====================================================== */

int moebius_mu (int n)
{
	if (1 >= n) return 1;
	if (3 >= n) return -1;

	INIT_PRIME_SIEVE(n);

	/* Implement the dumb/simple moebius algo */
	int cnt = 0;
	int i;
	for (i=0; ; i++)
	{
		int k = sieve[i];
		if (0 == n%k)
		{
			cnt ++;
			n /= k;

			/* If still divisible, its a square */
			if (0 == n%k) return 0;
		}

		if (1 == n) break;
		if (k*k > n)
		{
			cnt ++;
			break;
		}
	}

	if (0 == cnt%2) return 1;
	return -1;
}

/* ====================================================== */

int mertens_m (int n)
{
	int i;
	int acc = 0;
	for (i=1; i<=n; i++)
	{
		acc += moebius_mu (i);
	}
	return acc;
}

/* ====================================================== */

int carmichael_lambda (int n)
{
	if (9 >= n)
	{
		if (2 >= n) return 1;
		if (4 >= n) return 2;
		if (5 == n) return 4;
		if (6 == n) return 2;
		if (7 == n) return 6;
		if (8 == n) return 2;
		if (9 == n) return 6;
	}

	INIT_PRIME_SIEVE(n);

	/* Implement the dumb/simple factorization algo */
	int i;
	for (i=0; ; i++)
	{
		int p = sieve[i];
		if (0 == n%p)
		{
			int m = n;
			m /= p;
			while (0 == m%p) m /= p;

			// Case of Euler's function for prime powers.
			if (1 == m)
			{
				if (2 == p) return n/4; // weirdo special case.
				return (p-1)*n/p;
			}

			n = n / m;
			int t = (p-1)*n/p;
			if (2 == p) t = n/4;  // the special case, again.
			if (8 >= n) t = carmichael_lambda(n); // weird special case

			int r = carmichael_lambda(m);
			return lcm32(t,r);
		}
		// If we are here, and p*p > n,  then no prime less than
		// the sqrt of n divides it. Ergo, n must be prime.
		if (p*p > n) return n-1;
	}

	return 0;
}

/* ====================================================== */

int exp_mangoldt_lambda (int n)
{
	if (1 >= n) return 1;

	INIT_PRIME_SIEVE(n);

	/* Implement the dumb/simple factorization algo */
	int i;
	for (i=0; ; i++)
	{
		int p = sieve[i];
		if (0 == n%p)
		{
			n /= p;
			while (0 == n%p) n /= p;

			if (1 == n) return p;
			return 1;
		}
		// If we are here, and p*p > n,  then no prime less than
		// the sqrt of n divides it. Ergo, n must be prime.
		if (p*p > n) return n;
	}

	return 1;
}

/* ====================================================== */

long double mangoldt_lambda (int n)
{
	if (1 >= n) return 0.0L;

	int eml = exp_mangoldt_lambda(n);
	if (1 == eml) return 0.0L;
	return logl ((long double) eml);
}

/* ====================================================== */

long double mangoldt_lambda_cached (int n)
{
	DECLARE_LD_CACHE (mangoldt_cache);
	if (ld_one_d_cache_check (&mangoldt_cache, n))
	{
		return ld_one_d_cache_fetch(&mangoldt_cache, n);
	}
	else
	{
		long double val = mangoldt_lambda(n);
		ld_one_d_cache_store (&mangoldt_cache, val, n);
		return val;
	}
}

/* ====================================================== */

DECLARE_LD_CACHE (mangoldt_idx_cache);
DECLARE_UI_CACHE (mangoldt_pow_cache);
static int man_last_val =1;
static int man_last_idx =0;

long double mangoldt_lambda_indexed (int n)
{
	if (ld_one_d_cache_check (&mangoldt_idx_cache, n))
	{
		return ld_one_d_cache_fetch(&mangoldt_idx_cache, n);
	}
	else
	{
		ui_one_d_cache_check (&mangoldt_pow_cache, n);
		while (1)
		{
			man_last_val++;
			long double val = mangoldt_lambda(man_last_val);
			if (val != 0.0L)
			{
				man_last_idx++;
				ld_one_d_cache_store (&mangoldt_idx_cache, val, man_last_idx);
				ui_one_d_cache_store (&mangoldt_pow_cache, man_last_val, man_last_idx);
				if (n == man_last_idx)
					return val;
			}
		}
	}
}

unsigned int mangoldt_lambda_index_point (int n)
{
	if(ui_one_d_cache_check (&mangoldt_pow_cache, n))
	{
		return ui_one_d_cache_fetch(&mangoldt_pow_cache, n);
	}
	else
	{
		ld_one_d_cache_check (&mangoldt_idx_cache, n);
		while (1)
		{
			man_last_val++;
			long double val = mangoldt_lambda(man_last_val);
			if (val != 0.0L)
			{
				man_last_idx++;
				ld_one_d_cache_store (&mangoldt_idx_cache, val, man_last_idx);
				ui_one_d_cache_store (&mangoldt_pow_cache, man_last_val, man_last_idx);
				if (n == man_last_idx)
					return man_last_val;
			}
		}
	}
}

/* ====================================================== */
/* count the number of prime factors of n */

int liouville_omega (int n)
{
	int i;
	int acc = 2;

	if (1 >= n) return 0;
	if (2 >= n) return 1;

	INIT_PRIME_SIEVE(n);

	i=0;
	while (1)
	{
		int d = sieve[i];
		if (0 == n%d)
		{
			acc ++;
			n /= d;
		}
		else
		{
			i++;
		}
		if (d*d > n) break;
	}

	return acc-1;
}

int liouville_lambda (int n)
{
	int omega = liouville_omega (n);

	if (0 == omega%2) return 1;
	return -1;
}

/* ====================================================== */

#define TEST 1
#ifdef TEST

#include <stdio.h>

int test_isqrt (void)
{
	int have_error=0;
	int i;
	int nmax=10000;
	for (i=1; i<=nmax; i++)
	{
		int is = isqrt(i);
		if (is != integer_sqrt(i))
		{
			printf ("ERROR: in int sqrt function at n=%d\n", i);
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

/** Compute the divisor arithmetic function
 *  Returns the number of divisors of n.
 *  Raw brute force algorithm.
 */

int divisor_simple_algo (int n)
{
	int acc = 0;
	int d;

	for (d=1; d<= n; d++)
	{
		if (n%d) continue;
		acc ++;
	}

	return acc;
}

int test_divisor (void)
{
	int have_error=0;
	int i;
	int nmax=10000;
	for (i=1; i<=nmax; i++)
	{
		if (divisor(i) != divisor_simple_algo(i))
		{
			printf ("ERROR: in divisor function at n=%d\n", i);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested divisor function up to %d\n", nmax);
	}
	return have_error;
}

int test_moebius(void)
{
	int n;

	int have_error = 0;
	int nmax = 40000;
	for (n=1; n<nmax; n++)
	{
		/* Perform a Dirichlet sum */
		int sum = 0;
		int d;
		for (d=1; ; d++)
		{
			if (2*d > n) break;
			if (n%d) continue;
			sum += moebius_mu (d);
			// printf ("%d divides %d and sum=%d\n", d, n, sum);
		}
		if (1 != n) sum += moebius_mu (n);
		if (0 != sum)
		{
			printf ("ERROR for moebius mu at n=%d sum=%d \n", n, sum);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested moebius function w/ dirichlet up to %d\n", nmax);
	}
	return have_error;
}

int test_omega(void)
{
	int n;

	int have_error = 0;
	int nmax = 40000;
	for (n=1; n<nmax; n++)
	{
		/* Perform a Dirichlet sum */
		int sum = 0;
		int d;
		for (d=1; ; d++)
		{
			if (2*d > n) break;
			if (n%d) continue;
			sum += liouville_lambda (d);
			// printf ("%d divides %d and sum=%d\n", d, n, sum);
		}
		if (1 != n) sum += liouville_lambda (n);

		if (0 == sum) continue;

		int ns = sqrt (n);
		if (ns*ns != n)
		{
			printf ("ERROR at liouville lambda at n=%d sum=%d \n", n, sum);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested liouville omega function w/ dirichlet up to %d\n", nmax);
	}
	return have_error;
}

int main()
{
	test_isqrt();
	test_divisor ();
	test_omega ();
	test_moebius ();

	return 1;
}
#endif /* TEST */


/* --------------------------- END OF FILE ------------------------- */
