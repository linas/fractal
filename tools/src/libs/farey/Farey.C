
/* 
 * FUNCTION:
 * Convert real numbers to Farey numbers, and back.
 * Coded in psuedo-OO style for easy conversion to C++
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 * Linas Added stuff April 1996
 */

#ifdef LINUX
#include <signal.h>
#endif /* LINUX */

#include "Farey.h"
#include <stdio.h>
#include <math.h>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* ------------------------------------------------------------ */
ContinuedFraction::ContinuedFraction (void)
{
   real = 0.0;
   cutoff = 0x7fffffff;

   scratch = 0x0;

#ifdef LINUX
   /* block floating point exceptions */
   sigsetmask (0x180);
#endif /* LINUX */

   SetReal (0.0);
}

/* ------------------------------------------------------------ */

void 
ContinuedFraction::SetReal (double val) 
{
   RealToContinuedFraction (val);
}

/* ------------------------------------------------------------ */

void 
ContinuedFraction::SetRatio (int n, int d)
{
   RatioToContinuedFraction (n,d);
}

/* ------------------------------------------------------------ */

void 
ContinuedFraction::RatioToContinuedFraction (int numer, int deno)
{
   double tmp;
   int i;
   unsigned int n, d, m;

   real = ((double) numer) / ((double) deno);
   intpart = numer / deno;
   if (0.0 > real) intpart--;
   numer -= deno * intpart;
   num = numer;
   denom = deno;

   /* Next, compute the continued fraction of the ratio */
   for (i=0; i<32; i++) tinued_frac[i] = 0;
   n = num;
   d = denom;
   if (n != 0) {
      for (i=0; i<32; i++) {
         m = d/n;     
         tinued_frac [i] = m;
         /* printf (" term %d  denom = %d num = %d termval = %d ",i, d, n, m); */
         d -= m*n;
         /* printf ("rem = %d \n", d); */
         m = d;
         d = n;
         n = m;
   
         /* if (d/32>n)  */
         if (d/cutoff>n) {
            n = 0;   /* clamp for "virtually zero" */
            if (tinued_frac [i] == 1) {
               if (i != 0) {
                  tinued_frac [i] = 0;
                  tinued_frac [i-1] ++;
                  i--;
               } else {
                  intpart ++;
               }
            }
         }
   
         /* Check for termination of the expansion, and 
          * make sure that the continued fraction always has 
          * an even number of terms.  Huh?? Why do we do this?? */
         if (n == 0) {
            if (i%2 == 0) {
               tinued_frac [i] -= 1;
               tinued_frac [i+1] = 1;
            }
            break;
         }
      }
   }

   /* lets count the number of terms */
   for (i=0; i<32; i++) if (tinued_frac[i] == 0) break;
   nterms = i;
}

/* ------------------------------------------------------------ */

void 
ContinuedFraction::RealToContinuedFraction (double val)
{
   double tmp;
   unsigned int n, d, m;

   real = val;
   /* first, compute the integer and fractional parts of the real 
    * number x. Next, express the fractional part as a ratio of integers. */
   intpart = (int) val;
   if (0.0 > val) intpart--;
   fracpart = real - (double) intpart;
   d = 0x7fffffff;
   tmp = (double) d;
   tmp *= fracpart;
   n = (int) tmp;

   /* avoid rounding problems for floating point numbers */
   cutoff = 54321;
   RatioToContinuedFraction (n, d);
   cutoff = 0x7fffffff;

   intpart = (int) val;
   if (0.0 > val) intpart--;
}

/* ------------------------------------------------------------ */

void 
ContinuedFraction::Print (void)
{
   int i;
   printf (" ratio %d over %d is continued fraction of %d terms\n", num, denom, nterms);
   for (i=0; i<nterms; i++) { 
      printf (" term %d is %d partial=%g\n", i, tinued_frac[i], partial[i]);
   }
}

/* ------------------------------------------------------------ */

int 
ContinuedFraction::GetNumTerms (void)
{
   return nterms;
}

/* ------------------------------------------------------------ */
/* compute Farey Number from continued fraction.
 * Algorithm used here does base2 arithmetic by bit-shifting. 
 * For generalized algorithm, see below. 
 */

double 
ContinuedFraction::ToFarey (void)
{
   double tmp;
   int first_term;
   int i, j, k;
   unsigned int f, o;

   /* twidle the first bit */
   if (num != 0) {
      first_term = tinued_frac[0];
      tinued_frac[0] -= 1;
   }
   
   /* 
   printf ("\n\n"); 
   for (i=0; i<5; i++) { 
      printf (" term %d is %d \n", i, tinued_frac[i]);
   }
   printf ("\n\n");
   */

   /* now, build the binary rep of the Farey number */
   f = 0;
   o = 0x80000000;
   j = 0;
   for (i=0; i<32;) {
      if (j%2==0) {
         o >>= tinued_frac[j];
         /* printf (" tin = %d, o= 0x%x \n", tinued_frac[j], o); */
      } else {
         for (k=0; k<tinued_frac[j]; k++) {
            f |= o;
            o >>= 1;
         }
      }
      /* printf (" tin = %d f = 0x%x \n", tinued_frac[j], f); */
      i += tinued_frac[j];
      j++;
      if (tinued_frac[j] == 0) break;
   }

   /* finally, convert the farey number to floating point rep */
   tmp = 1.0 / (double) 0x7fffffff;
   tmp *= 0.5;   /* since 0xffffffff hangs up on minus sign */
   tmp *= (double) f;

   tmp += (double) intpart;

   /* reset the first term */
   tinued_frac[0] = first_term;

   return (tmp);
}

/* ------------------------------------------------------------ */
/* Converts continued fraction back to a real number */

double 
ContinuedFraction::ToReal (void)
{
   int i;
   double tmp;

   tmp = (double) intpart;
   real = tmp;
   if (nterms == 0) return (tmp);

   /* now, work backwards and reconstruct the fraction. */
   tmp = 1.0 / ((double) tinued_frac[nterms-1]);
   for (i=nterms-2; i>=0; i--) {
      tmp += (double) tinued_frac[i];
      tmp = 1.0 / tmp;
   }

   real += tmp;
   return (real);

}

/* ------------------------------------------------------------ */
/* Converts continued fraction into real number, 
 * but with a numerator of z instead of 1.
 */

double 
ContinuedFraction::ToZReal (double z)
{
   int i;
   double tmp;
   double znum;

   znum = (double) intpart;
   if (nterms == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   tmp = z / ((double) tinued_frac[nterms-1]);
   for (i=nterms-2; i>=0; i--) {
      tmp += (double) tinued_frac[i];
      /* not normally needed, will this help mystery crash? */
      if (tmp == 0.0) { tmp = 1.0e30; } else { tmp = z / tmp; }
   }

   /* get rid of last z, to normalize */
   tmp /= z;

   znum += tmp;
   return (znum);

}

/* ------------------------------------------------------------ */
/* Explicitly computes the size of the gaps in the z-expansion
 * basically, this is d/dz of ToZReal taken at z=1.0 
 */

double 
ContinuedFraction::ToZRealGap (void)
{
   int i;
   double tmp, deltahi, deltalo;

   if (nterms == 0) return (0);

   /* now, work backwards and reconstruct the fraction. */
   deltalo = 0.0;
   deltahi = 1.0;
   partial[nterms-1] = ((double) tinued_frac[nterms-1]);
   for (i=nterms-2; i>=0; i--) {
      tmp = 1.0 / partial[i+1];
      partial[i] = ((double) tinued_frac[i]) + tmp;
      deltahi = (1.0 - deltahi*tmp)*tmp;
      deltalo = (1.0 - deltalo*tmp)*tmp;
   }
   tmp = 1.0 / partial[0];
   deltahi = (1.0 - deltahi*tmp)*tmp;
   deltalo = (1.0 - deltalo*tmp)*tmp;

   // return deltalo;
   // return deltahi;
   return 0.5* (deltahi-deltalo);

}

/* ------------------------------------------------------------ */
/* converts continued fraction into real number, with a numerator of z.
 * inverts the whatever.
 */

double 
ContinuedFraction::ToInvZReal (double z)
{
   int i;
   double tmp;
   double znum;

   znum = (double) intpart;
   if (nterms == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   tmp = z * ((double) (tinued_frac[nterms-1]));
   for (i=nterms-2; i>=0; i--) {
      tmp += 1.0 / ((double) (tinued_frac[i]));
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

double 
ContinuedFraction::ToEReal (double t)
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) tinued_frac[i];
      tmp = exp (- ((double) i) * t) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an e-fraction to be the result of the continued fraction
 * where each term is damped by an exponential. 
 * Similar but slightly different than the above. Read the code. */

double 
ContinuedFraction::ToEFraction (double t)
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = exp (((double) (n-1)) * t);
   tmp *= ((double) tinued_frac[n-1]);
   tmp = 1.0 / tmp;
   for (i=n-2; i>=0; i--) {
      tmp += (double) tinued_frac[i] * exp (((double) i) * t);
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

double 
ContinuedFraction::ToEFarey (double t)
{
   int i;
   double sgn, sum, tmp;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<nterms; i++) {
      sum += (double) tinued_frac[i];
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
 * (Ooops. It not continuous at all rationals. It is continuous at
 * all rationals whose denominators are powers of two. Sigh. 
 */

double 
ContinuedFraction::ToSinFarey (double t)
{
   int i;
   double sgn, sum, tmp;
   double retval;
   double broke;
   double ct, cti;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<nterms; i++) {
      sum += (double) tinued_frac[i];
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

double 
ContinuedFraction::ToTFarey (double t)
{
   int i;
   double sgn, sum, tmp;
   double retval;
   double broke;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   tmp = 0.0;
   sum = -1.0;    /* initally, shift by one bit */
   sgn = 1.0;
   for (i=0; i<nterms; i++) {
      sum += (double) tinued_frac[i];
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

double 
ContinuedFraction::ToEFareyGap (double t)
{
   int i;
   double tmp, sum;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   sum = -2.0;    /* initally, shift by one bit */
   for (i=0; i<nterms; i++) {
      sum += (double) tinued_frac[i];
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

double 
ContinuedFraction::ToCosReal (double omega)
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) tinued_frac[i];
      tmp = cos (((double) i) * omega) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an bessel-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

double 
ContinuedFraction::ToSincReal (double omega)
{
   int i, n;
   double tmp;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) tinued_frac[i];
      tmp = j0 (((double) i) * omega) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

double 
ContinuedFraction::ToCnReal (double omega)
{
   int i, n;
   double tmp;
   double fomega;
   double retval;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) tinued_frac[i];
/*
      fomega = ((double) i) * omega;
      fomega *= 2.0 * M_PI;
      fomega = 2.0 * M_PI * omega / tmp; 
*/
      fomega = 2.0 * M_PI * omega;
      tmp = (1.0 + cos (fomega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

double 
ContinuedFraction::ToSnReal (double omega)
{
   int i, n;
   double tmp;
   double retval;
   double fomega;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {
      tmp += (double) tinued_frac[i];
/*
      fomega = ((double) i) * omega;
      fomega *= 2.0 * M_PI;
      fomega = 2.0 * M_PI * omega / tmp; 
*/
      fomega = 2.0 * M_PI * omega;
      tmp = (1.0 + sin (fomega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

double 
ContinuedFraction::ToFCnReal (double omega)
{
   int i, j, k, n;
   double tmp;
   double retval;
   double fomega;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   if (!scratch) scratch = new ContinuedFraction;

omega *= 2.0 * M_PI;
   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {

/*
      k = i;
*/
      k = 1;
      for (j=0; j<i; j++) k *= 3;

      scratch->SetReal (((double) k) * omega);
      fomega = scratch->ToFarey ();
/*
      fomega *= 2.0 * M_PI;
*/

      tmp += (double) tinued_frac[i];
      tmp = (1.0 + cos (fomega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an cos-real to be the result of the continued fraction
 * where each term is damped by a cosine. Read the code. */

double 
ContinuedFraction::ToFSnReal (double omega)
{
   int i, j, k, n;
   double tmp;
   double retval;
   double fomega;

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   if (!scratch) scratch = new ContinuedFraction;

omega *= 2.0 * M_PI;
   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {

/*
      k = i;
*/
      k = 1;
      for (j=0; j<i; j++) k *= 3;

      scratch->SetReal (((double) k) * omega);
      fomega = scratch->ToFarey();
/*
      fomega *= 2.0 * M_PI;
*/

      tmp += (double) tinued_frac[i];
      tmp = (1.0 + sin (fomega)) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */

double 
ContinuedFraction::ToZCnReal (double omega, double zzz)
{
   int i, j, k, n;
   double tmp;
   double retval;
   double fomega;
   double zees[101];

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   n = nterms;
   if (100 < n) n = 100;
   for (i=0; i<n; i++) {
      fomega = InvZReal (((double) i) * omega, zzz);
      fomega *= 2.0 * M_PI;
      zees[i] = fomega;
   }

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {

      tmp += (double) tinued_frac[i];
      tmp = (1.0 + cos (zees[i])) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */

double 
ContinuedFraction::ToZSnReal (double omega, double zzz)
{
   int i, j, k, n;
   double tmp;
   double retval;
   double fomega;
   double zees[101];

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   n = nterms;
   if (100 < n) n = 100;
   for (i=0; i<n; i++) {
      fomega = InvZReal (((double) i) * omega, zzz);
      fomega *= 2.0 * M_PI;
      zees[i] = fomega;
   }

   /* now, work backwards and reconstruct the fraction. */
   n = nterms;
   tmp = 0.0;
   for (i=n-1; i>=0; i--) {

      tmp += (double) tinued_frac[i];
      tmp = (1.0 + sin (zees[i])) / tmp;
   }

   retval += tmp;
   return (retval);
}

/* ------------------------------------------------------------ */
/* I've defined an x-fraction to be the result of the continued fraction
 * where each term is whatever. Read the code. */

double 
ContinuedFraction::ToXFarey (double t)
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

   retval = (double) intpart;
   if (nterms == 0) return (retval);

   logotoo = - log (2.0);
   sum = -1.0;
   cygnus = 1.0;

   /* now, work forwards */
   n = nterms;
   for (i=0; i< n; i++) {
      sum += ((double) tinued_frac[i]);
#ifdef BAD_STRAT
      if (i+1<n) {
         tmp = ((double) tinued_frac[i+1]);
         if (tmp > 1.1) {
            tmp *= t * exp (tmp * logotoo);
            tmp = 1.0 - tmp;
         } else{
            tmp = 1.0;
         }
      } else { tmp = 1.0;}
#endif BAD_STRAT
      if (0 < i) {
         tmp = ((double) tinued_frac[i-1]);
         if (tmp < 1.1) {
            tmp *= t * exp (((double) tinued_frac[i]) * logotoo);
            tmp = 1.0 - tmp;
         } else{
            tmp = 1.0;
         }
      } else { tmp = 1.0;}
      retval += cygnus * tmp * exp (sum * logotoo);
      cygnus *= -1.0;
/*
printf (" its %i %i %f  %f %f \n", i, tinued_frac[i], sum, tmp, retval);
*/
   }

   return (retval);
}

/* ------------------------------------------------------------ */
/* computes crazy sum of continued fractions
 */


double 
ContinuedFraction::CFSum (ContinuedFraction *other,
              double alpha, double beta, double gamma)
{
   int i;
   double eval;
   double znum;
   int minterms, maxterms;
   ContinuedFraction *tmp;

   znum = alpha * (double) intpart;
   znum += beta * ((double) (other->intpart));
   znum += gamma;

   minterms = MIN ((nterms), (other->nterms));
   maxterms = MAX ((nterms), (other->nterms));

   if (maxterms == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   eval = 0.0;
   if (maxterms == nterms) {
      for (i=nterms-1; i>=minterms; i--) {
         eval += alpha * ((double) tinued_frac[i]) + gamma;
         /* not normally needed, will this help mystery crash? */
         if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
      }
   } else {
      for (i=nterms-1; i>=minterms; i--) {
         eval += beta * ((double) (other->tinued_frac[i])) + gamma;
         /* not normally needed, will this help mystery crash? */
         if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
      }
   }

   for (i=minterms-1; i>=0; i--) {
      eval += alpha * ((double) tinued_frac[i]) + beta * ((double) other->tinued_frac[i]) + gamma;
      /* not normally needed, will this help mystery crash? */
      if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
   }

   znum += eval;
   return (znum);

}

/* ------------------------------------------------------------ */
/* computes crazy product of continued fractions
 */

double 
ContinuedFraction::CFProd (ContinuedFraction *other,
              double alpha, double beta)
{
   int i;
   double eval;
   double znum;
   int minterms, maxterms;
   struct Farey *tmp;

   znum = alpha * (double) intpart;
   znum += beta * ((double) (other->intpart));

   minterms = MIN ((nterms), (other->nterms));

   if (minterms == 0) return (znum);

   /* now, work backwards and reconstruct the fraction. */
   eval = 0.0;
   for (i=minterms-1; i>=0; i--) {
      eval += alpha * ((double) tinued_frac[i]) * ((double) other->tinued_frac[i]) + beta;
      /* not normally needed, will this help mystery crash? */
      if (eval == 0.0) { eval = 1.0e30; } else { eval = 1.0 / eval; }
   }

   znum += eval;
   return (znum);

}

/* ------------------------------------------------------------ */
/* computes inverses by binary subdivision.
 */

double Inverse (void *cxt, double (*func)(void *, double), double val)
{
   double mid, guess;
   int ires;
   double delta;

   delta = 0.5;
   guess = 0.5;

   mid = (*func)(cxt, guess);

   for (ires=0; ires <55; ires ++) {
      delta *= 0.5;

      if (val > mid) {
         guess += delta;
      } else {
         guess -= delta;
      }
      mid = (*func) (cxt, guess);
   }
   
   return guess;
}

/* ------------------------------------------------------------ */
/* computes inverse of Zreal mapping */

struct InvZReal_s {
   ContinuedFraction fp;
   double zee;
};

double InvZReal_f (void * stru, double val) 
{
   struct InvZReal_s *sp;
   double retval;

   sp = (struct InvZReal_s *) stru;
    
   sp->fp.SetReal (val);
   // retval = sp->fp.ToZReal (sp->zee);
   retval = sp->fp.ToFarey ();

   return retval;
}

double InvZReal (double val, double zzz) 
{
   struct InvZReal_s s;
   double retval;
   int intpart;

   s.zee = zzz;

   intpart = (int) val;
   if (0.0 > val) intpart --;
   val -= (double) intpart;

   retval = Inverse (((void *) &s), InvZReal_f, val);

   retval += (double) intpart;

   return retval;
}

/* ---------------------- END OF FILE ------------------------- */
