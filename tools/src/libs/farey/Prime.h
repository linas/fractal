
/*
 * FUNCTION:
 * Number Theory Stuff
 *
 * HISTORY:
 * Linas Vepstas 23 January 1994
 */

#ifndef __PRIME_H__
#define __PRIME_H__

struct Prime
{
   int nprimes;          /* number of primes in table */
   int *primes;          /* pointer to table  */

};

extern struct Prime * CreatePrime (void);
extern void PrintPrime (struct Prime *self);
extern int GetPrime (struct Prime *self, int i);

#endif /* __PRIME_H__ */
