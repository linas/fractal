/* totient.h: Calculates the Euler totient function, phi. 
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

#ifdef   __cplusplus
#extern "C" {
#endif

int totient_phi(int n);

#ifdef   __cplusplus
};
#endif
