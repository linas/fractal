
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include "Farey.h"
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------ */
/* returns the "farey-number" */

double f_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;

   RatioToContinuedFraction (f, nume, deno);
   r = ContinuedFractionToEFarey (f, z);

   return (r);
}

/* ------------------------------------------------------------ */
/* returns the "z-number" */

double r_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;

   RatioToContinuedFraction (f, nume, deno);
   r = ContinuedFractionToZReal (f, z);

   return (r);
}

/* ------------------------------------------------------------ */
/* returns a symmetric version of the z-number.
 * (the gaps are symmetrically arranged)
 */

double sym_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;

   r = r_sub_x (f, nume, deno, z);

   r -= r_sub_x (f, deno-nume, deno, z);
   r +=1.0;
   r *= 0.5;

   return (r);
}

/* ------------------------------------------------------------ */
/* Returns a symmetric z-number, 
 * where the gap at x=0.5 is directly given by the argument z
 */

double sym_gap_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;

   r = sym_sub_x (f, nume, deno, (1.0 - 2.0*z) /(1.0 + 2.0*z));

   return (r);
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   struct Farey *f;
   double r, x, y, z, t;
   double v, vprev, gap;
   int i, n;
   int nume, deno;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();

   n = atoi (argv[1]);
   z = atof (argv[2]);

   r = 0.0;
   v = 0.0;
   gap = 0.0;
   deno = n* 6899 +1;
   for (i=0; i<n; i++){
      nume = 6899*i;

      x = ((double) (nume))/ ((double) deno);
/*

      t = - log (0.5 * (1.0+z));
      r = f_sub_x (f, nume, deno, t);

      t = - log (0.5 * (1.0-z));
      r += f_sub_x (f, nume, deno, t);

      r *= 0.5;
*/
   /* get the farey number */
   vprev = v;
   RatioToContinuedFraction (f, nume, deno);
   v = ContinuedFractionToEFarey (f, z);
   if ((v-vprev) < 0.0) gap -= v-vprev;
   r = gap + v;

      printf ("i %g f %g \n", x, r);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
