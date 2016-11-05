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

long thue_morse(long n)
{
   if (0 == n) return 0;
   if (1 == n) return 1;
   if (0 == n%2) return thue_morse (n/2);
   return (1-thue_morse ((n-1)/2));
}

/* ====================================================== */
/**
 * Raise p to the n'th power.  Fast bit-shift algo.
 */

long ipow(long p, long n)
{
	long acc = 1;
	long psq = p;

	for (int i=0; i<sizeof(long); i++)
	{
		if (n & 0x1) acc *= psq;
		psq *= psq;
		n >>= 1;
		if (0 == n) return acc;
	}

	return acc;
}


/* ====================================================== */

static unsigned long *sieve = NULL;
static unsigned long sieve_size = 0;
static unsigned long sieve_max = 0;

#define INIT_PRIME_SIEVE(N) \
	if (!sieve || sieve[sieve_max]*sieve[sieve_max] <(N)) {\
		init_prime_sieve(N); \
	}

/* Initialize and fill in a prime-number sieve.
 * Handles primes up to 2^64
 *
 * XXX Ths should be converted to use the code in prime.h and prime.c,
 * but for now, its here, so as to maximize compiler optimization.
 */
static void
init_prime_sieve (long long prod)
{
	unsigned long n, j;
	unsigned long nstart;
	unsigned long pos;

	if (sieve)
	{
		long long ss = sieve[sieve_max];
		if (ss*ss > prod) return;
	}

	unsigned long max = 1000.0+sqrt (prod);

	if (!sieve)
	{
		sieve = (unsigned long *) malloc (8192*sizeof (unsigned long));
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
			long p = sieve[j];
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
					sieve = (unsigned long *)realloc (sieve, sieve_size * sizeof (unsigned long));
				}
				break;
			}
		}
	}
	sieve_max = pos-1;

#if 0
	for (j=0; j<pos; j++)
	{
		printf ("its %ld %ld\n", j, sieve[j]);
	}
#endif
}

/* ====================================================== */
/** Compute the divisor arithmetic function
 *  Returns the number of divisors of n.
 *  Almost a tail-recursive algorithm.
 */

static long divisor_helper (long n, long last_prime_checked)
{
	long ip;

	if (1==n) return 1;
	if (2==n) return 2;

	INIT_PRIME_SIEVE(n);

	for (ip=last_prime_checked+1; ; ip++)
	{
		long d = sieve[ip];

		if (d*d >n) break;
		if (n%d) continue;

		long acc = 2;
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

long divisor (long n)
{
	return divisor_helper (n, -1);
}

/* ====================================================== */
/** 
 * Sigma arithmetic series, equals divisor arith series for a=0
 * Computes the divisors of n, raises each to the a'th power, and
 * returns thier sum.
 *
 * Should be fast, works directly off the factorization.
 */
static long sigma_helper (long n, long a, long last_prime_checked)
{
	long ip;

	if (1==n) return 1;
	if (2==n) return 1 + ipow(2, a);

	INIT_PRIME_SIEVE(n);

	for (ip=last_prime_checked+1; ; ip++)
	{
		long d = sieve[ip];

		if (d*d>n) break;
		if (n%d) continue;

		long acc = 2;
		n /=d;
		while (0 == n%d)
		{
			n /=d;
			acc ++;
		}
		long fac = (ipow(d, a*acc) -1) / (ipow(d, a) - 1);
		if (1==n) return fac;
		return fac * sigma_helper (n, a, ip);
	}

	return 1 + ipow(n, a);
}

long sigma (long n, long a)
{
	if (0 == a) return divisor_helper(n, -1);
	return sigma_helper (n, a, -1);
}

/* sigma, but for float-point power, and a slow algo. */
long double sigmaf (long n, long double a)
{
	long double acc = 0;
	long d;

	long ns = n / 2;
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;
		acc += powl(d, a);
	}
	acc += powl(n, a);

	return acc;
}

/* same as above, but includes a log */
long double sigmalog (long n, long double a)
{
	long double acc = 0;
	long d;

	long ns = n / 2;
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;
		acc += powl(d, a) * logl(d);
	}
	acc += powl(n, a) * logl(n);

	return acc;
}

/* Sloww brute force algo */
long sigma_unitary (long n, long k)
{
	long acc = 0;
	long d;

	long ns = n/2;
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;

		if (1 != gcf64(d, n/d)) continue;
		long dp = 1;
		for (long ia=0; ia<k; ia++) dp *= d;
		acc += dp;
	}

	long dp = 1;
	for (long ia=0; ia<k; ia++) dp *= d;
	acc += dp;
	return acc;
}

/* ====================================================== */
DECLARE_UL_CACHE (sigma_one_cache);

long sigma_one (long n)
{
	if (ul_one_d_cache_check (&sigma_one_cache, n))
		return ul_one_d_cache_fetch(&sigma_one_cache, n);

	long val = sigma(n, 1);
	ul_one_d_cache_store (&sigma_one_cache, val, n);
	return val;
}

DECLARE_UL_CACHE (partition_cache);

long partition (long n)
{
	if (0 == n) return 1;

	// This over-flows a 64-bit int at n=316
	if (316 < n) return 0;

	if (ul_one_d_cache_check (&partition_cache, n))
		return ul_one_d_cache_fetch(&partition_cache, n);

	long acc = 0;
	for (long k=0; k < n; k++)
	{
		acc += sigma_one(n-k) * partition(k);
	}
	acc /= n;

	ul_one_d_cache_store (&partition_cache, acc, n);
	return acc;
}

DECLARE_ULL_CACHE (partition_ll_cache);

unsigned __int128 partitionll (long n)
{
	if (0 == n) return 1;

	// This over-flows a 128-bit int at n=1249
	if (1249 < n) return 0;

	if (ull_one_d_cache_check (&partition_ll_cache, n))
		return ull_one_d_cache_fetch(&partition_ll_cache, n);

	unsigned __int128 acc = 0;
	for (int k=0; k < n; k++)
	{
		unsigned __int128 sig = sigma_one(n-k);
		unsigned __int128 part = partitionll(k);
		acc += sig * part;
#if 0
if (403 < n) {
int sb = 0; unsigned __int128 sss = sig;
int pb = 0; unsigned __int128 ssp = part;
int ab = 0; unsigned __int128 ssa = acc;
while (0 < sss) { sb++; sss>>=1; }
while (0 < ssp) { pb++; ssp>>=1; }
while (0 < ssa) { ab++; ssa>>=1; }
unsigned long lo = acc & 0xffffffffffffffff;
unsigned long hi = acc >> 64;
printf("duuude n=%ld k=%d sigb=%d pbit=%d abit=%d acc=%lx %lx\n", n,k,sb, pb, ab, hi, lo);
}
#endif
	}
	acc /= (unsigned __int128) n;

	ull_one_d_cache_store (&partition_ll_cache, acc, n);
	return acc;
}

/* ====================================================== */

long moebius_mu (long n)
{
	if (1 >= n) return 1;
	if (3 >= n) return -1;

	INIT_PRIME_SIEVE(n);

	/* Implement the dumb/simple moebius algo */
	long cnt = 0;
	long i;
	for (i=0; ; i++)
	{
		long k = sieve[i];
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

long mertens_m (long n)
{
	long i;
	long acc = 0;
	for (i=1; i<=n; i++)
	{
		acc += moebius_mu (i);
	}
	return acc;
}

/* ====================================================== */

long carmichael_lambda (long n)
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
	long i;
	for (i=0; ; i++)
	{
		long p = sieve[i];
		if (0 == n%p)
		{
			long m = n;
			m /= p;
			while (0 == m%p) m /= p;

			// Case of Euler's function for prime powers.
			if (1 == m)
			{
				if (2 == p) return n/4; // weirdo special case.
				return (p-1)*n/p;
			}

			n = n / m;
			long t = (p-1)*n/p;
			if (2 == p) t = n/4;  // the special case, again.
			if (8 >= n) t = carmichael_lambda(n); // weird special case

			long r = carmichael_lambda(m);
			return lcm64(t,r);
		}
		// If we are here, and p*p > n,  then no prime less than
		// the sqrt of n divides it. Ergo, n must be prime.
		if (p*p > n) return n-1;
	}

	return 0;
}

/* ====================================================== */

long exp_mangoldt_lambda (long n)
{
	if (1 >= n) return 1;

	INIT_PRIME_SIEVE(n);

	/* Implement the dumb/simple factorization algo */
	long i;
	for (i=0; ; i++)
	{
		long p = sieve[i];
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

long double mangoldt_lambda (long n)
{
	if (1 >= n) return 0.0L;

	long eml = exp_mangoldt_lambda(n);
	if (1 == eml) return 0.0L;
	return logl ((long double) eml);
}

/* ====================================================== */

long double mangoldt_lambda_cached (long n)
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
DECLARE_UL_CACHE (mangoldt_pow_cache);
static long man_last_val =1;
static long man_last_idx =0;

long double mangoldt_lambda_indexed (long n)
{
	if (ld_one_d_cache_check (&mangoldt_idx_cache, n))
	{
		return ld_one_d_cache_fetch(&mangoldt_idx_cache, n);
	}

	ul_one_d_cache_check (&mangoldt_pow_cache, n);
	while (1)
	{
		man_last_val++;
		long double val = mangoldt_lambda(man_last_val);
		if (val != 0.0L)
		{
			man_last_idx++;
			ld_one_d_cache_store (&mangoldt_idx_cache, val, man_last_idx);
			ul_one_d_cache_store (&mangoldt_pow_cache, man_last_val, man_last_idx);
			if (n == man_last_idx)
				return val;
		}
	}
}

unsigned long mangoldt_lambda_index_point (long n)
{
	if (ul_one_d_cache_check (&mangoldt_pow_cache, n))
	{
		return ul_one_d_cache_fetch(&mangoldt_pow_cache, n);
	}

	ld_one_d_cache_check (&mangoldt_idx_cache, n);
	while (1)
	{
		man_last_val++;
		long double val = mangoldt_lambda(man_last_val);
		if (val != 0.0L)
		{
			man_last_idx++;
			ld_one_d_cache_store (&mangoldt_idx_cache, val, man_last_idx);
			ul_one_d_cache_store (&mangoldt_pow_cache, man_last_val, man_last_idx);
			if (n == man_last_idx)
				return man_last_val;
		}
	}
}

/* ====================================================== */
/* count the number of distinct prime factors of n */

long little_omega (long n)
{
	if (1 >= n) return 0;
	if (2 >= n) return 1;

	INIT_PRIME_SIEVE(n);

	long i=0;
	long acc = 0;
	while (1)
	{
		long d = sieve[i];
		if (0 == n%d) acc ++;
		while (0 == n%d) n /= d;
		i++;
		if (d*d > n) break;
	}
	if (1 != n) acc++;

	return acc;
}

/* ====================================================== */
/* count the number of prime factors of n */

long big_omega (long n)
{
	if (1 >= n) return 0;
	if (2 >= n) return 1;

	INIT_PRIME_SIEVE(n);

	long acc = 0;
	long i=0;
	while (1)
	{
		long d = sieve[i];
		while (0 == n%d)
		{
			acc ++;
			n /= d;
		}
		i++;
		if (d*d > n) break;
	}
	if (1 != n) acc++;

	return acc;
}

long liouville_lambda (long n)
{
	long omega = big_omega (n);

	if (0 == omega%2) return 1;
	return -1;
}

/* ====================================================== */

// #define TEST 1
#ifdef TEST

#include <stdio.h>

long test_pow (void)
{
	long have_error=0;
	if (8 != ipow(2,3)) have_error++;
	if (9 != ipow(3,2)) have_error++;
	if (1024 != ipow(2,10)) have_error++;
	if (27 != ipow(3,3)) have_error++;
	if (81*16 != ipow(6,4)) have_error++;

	if (0 == have_error)
	{
		printf ("PASS: tested power function\n");
	}
	return have_error;
}

/**
 * Compute the divisor arithmetic function
 * Returns the number of divisors of n.
 * Raw brute force algorithm.
 */
long divisor_simple_algo (long n)
{
	long acc = 0;
	for (long d=1; d<= n; d++)
	{
		if (n%d) continue;
		acc ++;
	}

	return acc;
}

long test_divisor (void)
{
	long have_error=0;
	long nmax=10000;
	for (long i=1; i<=nmax; i++)
	{
		if (divisor(i) != divisor_simple_algo(i))
		{
			printf ("ERROR: in divisor function at n=%ld\n", i);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested divisor function up to %ld\n", nmax);
	}
	return have_error;
}

long test_sigma_zero (void)
{
	long have_error=0;
	long nmax=10000;
	for (long i=1; i<=nmax; i++)
	{
		if (sigma(i, 0) != divisor_simple_algo(i))
		{
			printf ("ERROR: in sigma-zero function at n=%ld\n", i);
			printf ("wanted %ld got %ld\n", divisor_simple_algo(i), sigma(i, 0));
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested sigma-zero function up to %ld\n", nmax);
	}
	return have_error;
}

/**
 * Compute the sigma arithmetic function
 * Returns the sum of powers of divisors of n.
 * Raw brute force algorithm. Slow.
 */
long sigma_simple_algo (long n, long a)
{
	long acc = 0;
	long d;

	long ns = n/2;
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;

		long dp = 1;
		for (long ia=0; ia<a; ia++) dp *= d;
		acc += dp;
	}

	long dp = 1;
	for (long ia=0; ia<a; ia++) dp *= n;
	acc += dp;

	return acc;
}

long test_sigma (void)
{
	long have_error=0;

	// 1681 = 41^4 and so sigma(1681,4) = (41^20 -1) / (41^4 -1) = 6e25
	// which overflows 64-bits .. almost overflows 128 bits.
	long nmax=1680;
	long smax = 5;
	for (long i=1; i<=nmax; i++)
	{
		for (long s=1; s<smax; s++)
		{
			if (sigma(i, s) != sigma_simple_algo(i, s))
			{
				printf ("ERROR: in sigma %ld function at n=%ld\n", s, i);
				printf ("wanted %ld got %ld\n", sigma_simple_algo(i, s), sigma(i, s));
				have_error ++;
			}
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested sigma function up a=%ld and n to %ld\n", smax, nmax);
	}
	return have_error;
}

long test_unitary_divisor (void)
{
	long have_error=0;
	long nmax=10000;
	for (long i=1; i<=nmax; i++)
	{
		long ud = 1 << little_omega(i);
		if (sigma_unitary(i, 0) != ud)
		{
			printf ("ERROR: in unitary divisor function at n=%ld\n", i);
			printf ("wanted %ld got %ld\n", ud, sigma_unitary(i, 0));
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested unitary divisor function up to %ld\n", nmax);
	}
	return have_error;
}

long test_partition (void)
{
	long have_error=0;

	if (101 != partition(13)) have_error++;
	if (77 != partition(12)) have_error++;
	if (56 != partition(11)) have_error++;
	if (42 != partition(10)) have_error++;
	if (30 != partition(9)) have_error++;
	if (22 != partition(8)) have_error++;
	if (15 != partition(7)) have_error++;
	if (11 != partition(6)) have_error++;
	if (7 != partition(5)) have_error++;
	if (5 != partition(4)) have_error++;
	if (3 != partition(3)) have_error++;
	if (2 != partition(2)) have_error++;
	if (1 != partition(1)) have_error++;
	if (1 != partition(0)) have_error++;

	long nmax=316;
	for (long i=101; i<=nmax; i++)
	{
		double asymp = exp(M_PI * sqrt(2.0*i/3.0)) / (4.0*i*sqrt(3));
		long ub = (long) asymp;
		long lb = (long) (0.95 * asymp);
		long part = partition(i);

		if (part < lb || ub < part)
		{
			int b = 0; long ss = part;
			while (0 < ss) { b++; ss>>=1; }
			printf ("ERROR: in paritition function at n=%ld\n", i);
			printf ("wanted %ld < %ld < %ld at bits=%d\n", lb, part, ub, b);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested parition function up to %ld\n", nmax);
	}
	else
	{
		printf ("FAIL: parition function is bad\n");
	}
	return have_error;
}

long test_partitionll (void)
{
	long have_error=0;

	if (101 != partitionll(13)) have_error++;
	if (77 != partitionll(12)) have_error++;
	if (56 != partitionll(11)) have_error++;
	if (42 != partitionll(10)) have_error++;
	if (30 != partitionll(9)) have_error++;
	if (22 != partitionll(8)) have_error++;
	if (15 != partitionll(7)) have_error++;
	if (11 != partitionll(6)) have_error++;
	if (7 != partitionll(5)) have_error++;
	if (5 != partitionll(4)) have_error++;
	if (3 != partitionll(3)) have_error++;
	if (2 != partitionll(2)) have_error++;
	if (1 != partitionll(1)) have_error++;
	if (1 != partitionll(0)) have_error++;

	long nmax=1249;
	for (long i=101; i<=nmax; i++)
	{
		double asymp = exp(M_PI * sqrt(2.0*i/3.0)) / (4.0*i*sqrt(3));
		unsigned __int128 ub = (unsigned __int128) asymp;
		unsigned __int128 lb = (unsigned __int128) (0.95 * asymp);
		if (400< i) lb = (unsigned __int128) (0.97 * asymp);
		if (800< i) lb = (unsigned __int128) (0.98 * asymp);
		unsigned __int128 part = partitionll(i);

		if (part < lb || ub < part)
		{
			int b = 0; unsigned __int128 ss = part;
			while (0 < ss) { b++; ss>>=1; }
			printf ("ERROR: in parititionll functionll at n=%ld\n", i);
			double flb = lb;
			double fub = ub;
			double fp = part;
			printf ("wanted %g < %g < %g at bits=%d\n", flb, fp, fub, b);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested paritionll function up to %ld\n", nmax);
	}
	else
	{
		printf ("FAIL: paritionll function is bad\n");
	}
	return have_error;
}

long test_moebius(void)
{
	long have_error = 0;
	long nmax = 40000;
	for (long n=1; n<nmax; n++)
	{
		/* Perform a Dirichlet sum */
		long sum = 0;
		long d;
		for (d=1; ; d++)
		{
			if (2*d > n) break;
			if (n%d) continue;
			sum += moebius_mu (d);
			// printf ("%d divides %ld and sum=%ld\n", d, n, sum);
		}
		if (1 != n) sum += moebius_mu (n);
		if (0 != sum)
		{
			printf ("ERROR for moebius mu at n=%ld sum=%ld \n", n, sum);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested moebius function w/ dirichlet up to %ld\n", nmax);
	}
	return have_error;
}

long test_lambda(void)
{
	long have_error = 0;
	long nmax = 40000;
	for (long n=1; n<nmax; n++)
	{
		/* Perform a Dirichlet sum */
		long sum = 0;
		for (long d=1; ; d++)
		{
			if (2*d > n) break;
			if (n%d) continue;
			sum += liouville_lambda (d);
			// printf ("%d divides %ld and sum=%ld\n", d, n, sum);
		}
		if (1 != n) sum += liouville_lambda (n);

		if (0 == sum) continue;

		long ns = sqrt (n);
		if (ns*ns != n)
		{
			printf ("ERROR at liouville lambda at n=%ld sum=%ld \n", n, sum);
			have_error ++;
		}
	}
	if (0 == have_error)
	{
		printf ("PASS: tested liouville lambda function w/ dirichlet up to %ld\n", nmax);
	}
	return have_error;
}

int main()
{
	test_pow ();
	test_divisor ();
	test_sigma_zero ();
	test_sigma ();
	test_unitary_divisor();
	test_partition ();
	test_partitionll ();
	test_lambda ();
	test_moebius ();

	return 1;
}
#endif /* TEST */


/* --------------------------- END OF FILE ------------------------- */
