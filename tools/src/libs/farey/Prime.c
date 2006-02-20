
/*
 * FUNCTION:
 * Number Theory Stuff
 *
 * HISTORY:
 * Linas Vepstas 23 January 1994
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Prime.h"


#define NPRIMES (self->nprimes)
#define PRIME (self->primes)

#define INIT_TABLE (1500)

/* ====================================================== */

struct Prime * CreatePrime () 
{
   struct Prime *self;
   int next;
   int i, j, max;

   self = (struct Prime *) malloc (sizeof (struct Prime *));

   NPRIMES = INIT_TABLE;
   PRIME = (int *)malloc (NPRIMES * sizeof (int));

   /* fill in the first few by hand */
   PRIME [0] = 2;
   PRIME [1] = 3;
   PRIME [2] = 5;
   PRIME [3] = 7;
   PRIME [4] = 11;
   PRIME [5] = 13;
   PRIME [6] = 17;
   PRIME [7] = 19;
   PRIME [8] = 23;
   PRIME [9] = 29;
   PRIME [10] = 31;

   next = 11;
   for (i=37; next < NPRIMES; i+=2) {
      
      /* first few by hand */
      if (i%3 == 0) continue;
      if (i%5 == 0) continue;
      if (i%7 == 0) continue;
      if (i%11 == 0) continue;
      if (i%13 == 0) continue;
      if (i%17 == 0) continue;
      if (i%19 == 0) continue;
      if (i%23 == 0) continue;
      if (i%29 == 0) continue;
      if (i%31 == 0) continue;

      /* the rest automatically */
      max = (int) sqrt ((double) i) +1;
      for (j=11; j<max; j++) {
         if (i%PRIME[j] == 0) continue;
      }

      /* if (next %100 ==0) { printf ("."); fflush (stdout); } */
      /* if we made it to here, it must be prime */
      PRIME[next] = i;
      next ++;
   }

   return (self);
}

/* ====================================================== */

void PrintPrime (struct Prime *self)
{
   int i;

   for (i=0; i< NPRIMES; i++) {
      printf ("The %d th prime is %d \n", i, PRIME[i]);
   }
   printf ("\n\n");
}

/* ====================================================== */

int GetPrime (struct Prime *self, int i)
{
   if (i < NPRIMES) {
      return (PRIME[i]);
   } else {
      return 0;
   }
}

/* ====================================================== */
/* Implements erasthmos seive test for primeness */

int IsPrime (struct Prime *self, int num) 
{
   double tmp;
   int max;

   tmp = num;
   tmp = sqrt (num);
   max = tmp;
   max ++;

	return 1;
}

/* ====================================================== */
void
test_prime ()
{
   struct Prime *p;

   p = CreatePrime ();

   PrintPrime (p);
}

/* ====================================================== */
