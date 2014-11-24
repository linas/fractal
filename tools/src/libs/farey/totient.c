/* totient.c: Calculates the Euler totient function, phi. 
 *
 *  By Terry R. McConnell
 *
 *                              Theory
 *
 *  Euler's totient function, phi(n), is defined as the number of natural
 *  numbers <= n which are relatively prime to n. (We count 1 as relatively
 *  prime to everything.) It is an extremely important function in number
 *  theory. For primes p it is clear from the definition that phi(p) = p - 1.
 *  For powers of a prime it also easy to see (use induction on n) that
 *  phi(p^n) = p^(n-1)(p-1). Thus, e.g, phi(125) = 100. For all other
 *  values phi can be computed by factoring n completely and using the
 *  following result:
 *
 *  Theorem: Phi is a multiplicative function, i.e., if (m,n) = 1 (relatively
 *  prime) then phi(mn) = phi(m)phi(n).
 *
 *  Perhaps the easiest proof is via group theory. It follows from the
 *  Chinese Remainder Theorem that the multiplicative group of invertible 
 *  integers modulo  mn, Z(mn), is isomorphic to the direct product Z(m)xZ(n).
 *  The result follows since the orders of these groups are given by the 
 *  totient function. 
 *
 *  For a proof of the version of the Chinese Remainder Theorem required,
 *  see Theorem 3.7 in H.M. Stark, An Introduction to Number Theory, 
 *  Markham, Chicago, 1970.
 *
 *  Our implementation is highly recursive and is motivated by the McCarthy
 *  conditional statement formalism. (See Marvin Minsky, Computation:
 *  finite and infinite machines, Prentice Hall, Englewood Cliffs, 1967,
 *  problem 10.7-1.)
 *
 *  We define phi(-n) = phi(n), and do not define phi(0).
 *
*/

#include "gcf.h"

/* In conditional statement form, phi can be defined together with another
 * function of 2 variables we denote as phiphi. We have phi(n) = phiphi(n,2)
 * and phiphi(y,x) = if(x = y-1 then x else if( x|y then
 * 	if((x,y/x)=1 then phi(x)phi(y/x) else x phi(y/x)) else phi(y,x+1))))
 */

static int phiphi(int,int);
int totient_phi(int n)
{
	if (n < 0) n=-n;
	/* handle a few trivial boundary cases */
	if (n <= 1) return 0;
	if (n == 2) return 1;
	if (n == 3) return 2;
	return phiphi(n,2);
}

/* This only gets called with y >= 3 and y > x >= 2 */

static int phiphi(int y, int x)
{
	int z;

	if (x+1 == y) return x; /* phi(prime p) = p-1 */
	if ((y%x) == 0)
	{
		z = y/x;
		if (gcf32(x,z) == 1)
			return totient_phi(x)*totient_phi(z); /* multiplicative property */
		else
			return x*totient_phi(z); /* This is a tricky case. It may
					    happen when x is a prime such
					    that a power of x divides y. In
					    case y = p^n, phi(y) = p^(n-1)(p-1)
					   */
	}
	else return phiphi(y,x+1);
}

