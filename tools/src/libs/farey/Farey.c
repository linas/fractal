
/* 
 * FUNCTION:
 * Convert real numbers to Farey numbers, and back
 * Coded in psuedo-OO style for easy conversion to C++
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 */

#ifdef LINUX
#include <signal.h>
#endif /* LINUX */

#include "Farey.h"
#include <stdio.h>
#include <math.h>

#define REAL (self->real)
#define INTPART (self->intpart)
#define FRACPART (self->fracpart)
#define NUM (self->num)
#define DENOM (self->denom)
#define TINUED_FRAC (self->tinued_frac)
#define NTERMS (self->nterms)
#define CUTOFF (self->cutoff)

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* ------------------------------------------------------------ */
struct Farey * CreateFarey ()
{
   struct Farey * self;
   self = (struct Farey *) malloc (sizeof (struct Farey));
   CUTOFF = 0x7fffffff;

#ifdef LINUX
   /* block floating point exceptions */
   sigsetmask (0xffff);
#endif /* LINUX */
   return (self);
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
void SetReal (struct Farey *self, double val) 
#else
void SetReal (self, val) 
struct Farey *self;
double val;
#endif
{
   RealToContinuedFraction (self, val);
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
double GetReal (struct Farey *self)
#else
double GetReal (self) 
struct Farey *self;
#endif
{
   ContinuedFractionToReal (self);
   return (REAL);
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
double GetFarey (struct Farey *self)
#else
double GetFarey (self) 
struct Farey *self;
#endif
{
   return ContinuedFractionToFarey (self);
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
void RatioToContinuedFraction (struct Farey *self, int numer, int deno)
#else
void RatioToContinuedFraction (self, numer, deno)
struct Farey *self;
int numer, deno;
#endif
{
   double tmp;
   int i;
   unsigned int n, d, m;

   REAL = ((double) numer) / ((double) deno);
   INTPART = numer / deno;
   numer -= deno * INTPART;
   NUM = numer;
   DENOM = deno;

   /* Next, compute the continued fraction of the ratio */
   for (i=0; i<32; i++) TINUED_FRAC[i] = 0;
   n = NUM;
   d = DENOM;
   if (n != 0) {
      for (i=0; i<32; i++) {
         m = d/n;     
         TINUED_FRAC [i] = m;
         /* printf (" term %d  denom = %d num = %d termval = %d ",i, d, n, m); */
         d -= m*n;
         /* printf ("rem = %d \n", d); */
         m = d;
         d = n;
         n = m;
   
         /* if (d/32>n)  */
         if (d/CUTOFF>n) {
            n = 0;   /* clamp for "virtually zero" */
            if (TINUED_FRAC [i] == 1) {
               if (i != 0) {
                  TINUED_FRAC [i] = 0;
                  TINUED_FRAC [i-1] ++;
                  i--;
               } else {
                  INTPART ++;
               }
            }
         }
   
         /* Check for termination of the expansion, and 
          * make sure that the continued fraction always has 
          * an even number of terms */
         if (n == 0) {
            if (i%2 == 0) {
               TINUED_FRAC [i] -= 1;
               TINUED_FRAC [i+1] = 1;
            }
            break;
         }
      }
   }

   /* lets count the number of terms */
   for (i=0; i<32; i++) if (TINUED_FRAC[i] == 0) break;
   NTERMS = i;
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
void RealToContinuedFraction (struct Farey *self, double val)
#else
void RealToContinuedFraction (self, val)
struct Farey *self;
double val;
#endif
{
   double tmp;
   unsigned int n, d, m;

   REAL = val;
   /* first, compute the integer and fractional parts of the real 
    * number x. Next, express the fractional part as a ratio of integers. */
   INTPART = (int) REAL;
   FRACPART = REAL - (double) INTPART;
   d = 0x7fffffff;
   tmp = (double) d;
   tmp *= FRACPART;
   n = (int) tmp;

   /* avoid rounding problems for floating point numbers */
   CUTOFF = 54321;
   RatioToContinuedFraction (self, n, d);
   CUTOFF = 0x7fffffff;
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
void PrintContinuedFraction (struct Farey *self)
#else
void PrintContinuedFraction (self)
struct Farey *self;
#endif
{
   int i;
   printf (" ratio %d over %d is continued fraction of %d terms\n", NUM, DENOM, NTERMS);
   for (i=0; i<NTERMS; i++) { 
      printf (" term %d is %d \n", i, TINUED_FRAC[i]);
   }
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
double ContinuedFractionTerms (struct Farey *self)
#else
double ContinuedFractionTerms (self)
struct Farey *self;
#endif
{
   return ((double) NTERMS);
}

/* ------------------------------------------------------------ */
/* compute Farey Number from continued fraction.
 * Algorithm used here does base2 arithmetic by bit-shifting. 
 * For generalized algorithm, see below. 
 */

#ifdef ANSI_C
double ContinuedFractionToFarey (struct Farey *self)
#else
double ContinuedFractionToFarey (self)
struct Farey *self;
#endif
{
   double tmp;
   int first_term;
   int i, j, k;
   unsigned int f, o;

   /* twidle the first bit */
   if (NUM != 0) {
      first_term = TINUED_FRAC[0];
      TINUED_FRAC[0] -= 1;
   }
   
   /* 
   printf ("\n\n"); 
   for (i=0; i<5; i++) { 
      printf (" term %d is %d \n", i, TINUED_FRAC[i]);
   }
   printf ("\n\n");
   */

   /* now, build the binary rep of the Farey number */
   f = 0;
   o = 0x80000000;
   j = 0;
   for (i=0; i<32;) {
      if (j%2==0) {
         o >>= TINUED_FRAC[j];
         /* printf (" tin = %d, o= 0x%x \n", TINUED_FRAC[j], o); */
      } else {
         for (k=0; k<TINUED_FRAC[j]; k++) {
            f |= o;
            o >>= 1;
         }
      }
      /* printf (" tin = %d f = 0x%x \n", TINUED_FRAC[j], f); */
      i += TINUED_FRAC[j];
      j++;
      if (TINUED_FRAC[j] == 0) break;
   }

   /* finally, convert the farey number to floating point rep */
   tmp = 1.0 / (double) 0x7fffffff;
   tmp *= 0.5;   /* since 0xffffffff hangs up on minus sign */
   tmp *= (double) f;

   tmp += (double) INTPART;

   /* reset the first term */
   TINUED_FRAC[0] = first_term;

   return (tmp);
}

/* ------------------------------------------------------------ */

#ifdef ANSI_C
double ContinuedFractionToReal (struct Farey *self)
#else
double ContinuedFractionToReal (self)
struct Farey *self;
#endif
{
   int i;
   double tmp;

   tmp = (double) INTPART;
   REAL = tmp;
   if (NTERMS == 0) return (tmp);

   /* now, work backwards and reconstruct the fraction. */
   tmp = 1.0 / ((double) TINUED_FRAC[NTERMS-1]);
   for (i=NTERMS-2; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = 1.0 / tmp;
   }

   REAL += tmp;
   return (REAL);

}

/* ------------------------------------------------------------ */
/* converts continued fraction into real number, with a numerator of z.
 */


#ifdef ANSI_C
double ContinuedFractionToZReal (struct Farey *self, double z)
#else
double ContinuedFractionToZReal (self, z)
struct Farey *self;
double z;
#endif
{
   int i;
   double tmp;
   double znum;

   znum = (double) INTPART;
   if (NTERMS == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   tmp = z / ((double) TINUED_FRAC[NTERMS-1]);
   for (i=NTERMS-2; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      /* not normally needed, will this help mystery crash? */
      if (tmp == 0.0) { tmp = 1.0e30; } else { tmp = z / tmp; }
   }

   /* get rid of last z, to normalize */
   tmp /= z;

   znum += tmp;
   return (znum);

}

/* ------------------------------------------------------------ */
/* I've defined an e-real to be the result of the continued fraction
 * where each term is damped by an exponential. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToEReal (struct Farey *self, double t)
#else
double ContinuedFractionToEReal (self, t)
struct Farey *self;
double t;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = exp (- ((double) i) * t) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an e-fraction to be the result of the continued fraction
 * where each term is damped by an exponential. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToEFraction (struct Farey *self, double t)
#else
double ContinuedFractionToEFraction (self, t)
struct Farey *self;
double t;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = exp (((double) (n-1)) * t);
   tmp *= ((double) TINUED_FRAC[n-1]);
   tmp = 1.0 / tmp;
   for (i=n-2; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i] * exp (((double) i) * t);
      tmp = 1.0 / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* the routine below computes the generalized farey number from the
 * continued fraction. To get the usual farey number, simply use t=ln2;
 * this gives same result as bit-shifting algorithm above.
 */

#ifdef ANSI_C
double ContinuedFractionToEFarey (struct Farey *self, double t)
#else
double ContinuedFractionToEFarey (self, t)
struct Farey *self;
double t;
#endif
{
   int i;
   double sgn, sum, tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<NTERMS; i++) {
      sum += (double) TINUED_FRAC[i];
      tmp += sgn * exp (-t*sum);
      sgn = - sgn;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* the routine below computes the generalized farey number from the
 * continued fraction. To get the usual farey number, simply use t=1;
 * this gives same result as bit-shifting algorithm above.
 *
 * The algorithm multiplies each term by g sub k (n) where
 * g sub k (n) = (1 + sin ((k+1)*t) / sin (k*t)) ** (-n)
 * If I computed this correctly, this function is continous at all rationals.
 */

#ifdef ANSI_C
double ContinuedFractionToSinFarey (struct Farey *self, double t)
#else
double ContinuedFractionToSinFarey (self, t)
struct Farey *self;
double t;
#endif
{
   int i;
   double sgn, sum, tmp;
   double retval;
   double broke;
   double ct, cti;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<NTERMS; i++) {
      sum += (double) TINUED_FRAC[i];
      ct = 0.5 * (1.0 + cos (((double) i)*t));
      cti = 0.5 * (1.0 + cos (((double)(i+1))*t));
      broke = 1.0 + cti/ct;
      tmp += sgn * ct * pow (broke, -sum);
      sgn = - sgn;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* the routine below computes the generalized farey number from the
 * continued fraction. To get the usual farey number, simply use t=1;
 * this gives same result as bit-shifting algorithm above.
 *
 * The algorithm multiplies each term by g sub k (n) where
 * g sub k (n) = t**k / (1+t)**n
 * If I computed this correctly, this function is continous at all rationals.
 */

#ifdef ANSI_C
double ContinuedFractionToTFarey (struct Farey *self, double t)
#else
double ContinuedFractionToTFarey (self, t)
struct Farey *self;
double t;
#endif
{
   int i;
   double sgn, sum, tmp;
   double retval;
   double broke;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<NTERMS; i++) {
      sum += (double) TINUED_FRAC[i];
      broke = sum * log (1.0 + t);          /* linux pow() is broken */
      if (600.0 < broke ) broke = 600.0;
      tmp += sgn * pow (t, ((double) i)) * exp (-broke);
      sgn = - sgn;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* the routine below returns the "gap" between generalized Farey Numbers
 */

#ifdef ANSI_C
double ContinuedFractionToEFareyGap (struct Farey *self, double t)
#else
double ContinuedFractionToEFareyGap (self, t)
struct Farey *self;
double t;
#endif
{
   int i;
   double tmp, sum;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   sum = -2.0;    /* initally, shift by one bit */
   for (i=0; i<NTERMS; i++) {
      sum += (double) TINUED_FRAC[i];
   }

   if (sum > 0.0) {
      tmp = exp (-t * sum);
   } else {
      tmp = 1.0;
   }
   tmp *= 1.0 - 2.0 * exp (-t);
   retval += tmp;
   return (retval);
}


/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToCosReal (struct Farey *self, double omega)
#else
double ContinuedFractionToCosReal (self, omega)
struct Farey *self;
double omega;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = cos (((double) i) * omega) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an bessel-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToSincReal (struct Farey *self, double omega)
#else
double ContinuedFractionToSincReal (self, omega)
struct Farey *self;
double omega;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = j0 (((double) i) * omega) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToCnReal (struct Farey *self, double omega)
#else
double ContinuedFractionToCnReal (self, omega)
struct Farey *self;
double omega;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = (1.0 + cos (((double) i) * omega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToSnReal (struct Farey *self, double omega)
#else
double ContinuedFractionToSnReal (self, omega)
struct Farey *self;
double omega;
#endif
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = NTERMS;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) TINUED_FRAC[i];
      tmp = (1.0 + sin (((double) i) * omega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an x-fraction to be the result of the continued fraction
 * where each term is whaever. Read the code. */

#ifdef ANSI_C
double ContinuedFractionToXFarey (struct Farey *self, double t)
#else
double ContinuedFractionToXFarey (self, t)
struct Farey *self;
double t;
#endif
{
   int i, n;
   double tmp;
   double retval;
   double logotoo;
   double cygnus;
   double sum;
/*
printf ("\n\n yooo \n");
*/

   retval = (double) INTPART;
   if (NTERMS == 0) return (retval);

   logotoo = - log (2.0);
   sum = -1.0;
   cygnus = 1.0;

   /* now, work forwards */
   n = NTERMS;
   for (i=0; i< n; i++) {
      sum += ((double) TINUED_FRAC[i]);
#ifdef BAD_STRAT
      if (i+1<n) {
         tmp = ((double) TINUED_FRAC[i+1]);
         if (tmp > 1.1) {
            tmp *= t * exp (tmp * logotoo);
            tmp = 1.0 - tmp;
         } else{
            tmp = 1.0;
         }
      } else { tmp = 1.0;}
#endif BAD_STRAT
      if (0 < i) {
         tmp = ((double) TINUED_FRAC[i-1]);
         if (tmp < 1.1) {
            tmp *= t * exp (((double) TINUED_FRAC[i]) * logotoo);
            tmp = 1.0 - tmp;
         } else{
            tmp = 1.0;
         }
      } else { tmp = 1.0;}
      retval += cygnus * tmp * exp (sum * logotoo);
      cygnus *= -1.0;
/*
printf (" its %i %i %f  %f %f \n", i, TINUED_FRAC[i], sum, tmp, retval);
*/
   }

   return (retval);
}

/* ------------------------------------------------------------ */
/* computes crazy sum of continued fractions
 */


#ifdef ANSI_C
double CFSum (struct Farey *self, struct Farey *other,
              double alpha, double beta, double gamma)
#else
double CFSum (self, other, alpha, beta, gamma)
struct Farey *self, *other;
double alpha, beta, gamma;
#endif
{
   int i;
   double eval;
   double znum;
   int minterms, maxterms;
   struct Farey *tmp;

   znum = alpha * (double) INTPART;
   znum += beta * ((double) (other->intpart));
   znum += gamma;

   minterms = MIN ((NTERMS), (other->nterms));
   maxterms = MAX ((NTERMS), (other->nterms));

   if (maxterms == 0) return (znum);

   /* let self point at the one with more terms */
   /* this s NOT an object-oriented move. HACK_LAERT */
   if (maxterms != NTERMS) {
      tmp = other;
      other = self;
      self = tmp;

      eval = alpha;
      alpha = beta;
      beta = eval;
   }

   /* now, work backwards and reconstruct the fraction. */
   eval = 0.0;
   for (i=NTERMS-1; i>=minterms; i--) {
      eval += alpha * ((double) TINUED_FRAC[i]) + gamma;
      /* not normally needed, will this help mystery crash? */
      if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
   }

   for (i=minterms-1; i>=0; i--) {
      eval += alpha * ((double) TINUED_FRAC[i]) + beta * ((double) other->tinued_frac[i]) + gamma;
      /* not normally needed, will this help mystery crash? */
      if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
   }

   znum += eval;
   return (znum);

}

/* ------------------------------------------------------------ */
/* computes crazy product of continued fractions
 */


#ifdef ANSI_C
double CFProd (struct Farey *self, struct Farey *other,
              double alpha, double beta)
#else
double CFProd (self, other, alpha, beta)
struct Farey *self, *other;
double alpha, beta;
#endif
{
   int i;
   double eval;
   double znum;
   int minterms, maxterms;
   struct Farey *tmp;

   znum = alpha * (double) INTPART;
   znum += beta * ((double) (other->intpart));

   minterms = MIN ((NTERMS), (other->nterms));

   if (minterms == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   eval = 0.0;
   for (i=minterms-1; i>=0; i--) {
      eval += alpha * ((double) TINUED_FRAC[i]) * ((double) other->tinued_frac[i]) + beta;
      /* not normally needed, will this help mystery crash? */
      if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
   }

   znum += eval;
   return (znum);

}

/* ---------------------- END OF FILE ------------------------- */
