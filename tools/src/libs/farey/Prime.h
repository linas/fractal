
/*
 * FUNCTION:
 * Number Theory Stuff
 *
 * HISTORY:
 * Linas Vepstas 23 January 1994
 */

#include <math.h>

struct Prime {
   int nprimes;          /* number of primes in table */
   int *primes;          /* pointer to table  */

};

#ifdef ANSI_C

extern struct Prime * CreatePrime (void);
extern void PrintPrime (struct Prime *self);
extern int GetPrime (struct Prime *self, int i);

#else

extern struct Prime * CreatePrime ();
extern void PrintPrime ();
extern int GetPrime ();

#endif /* ANSI_C */

/* ====================================================== */
