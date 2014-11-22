/*
 * FUNCTION:
 * This module contains code for generating gaussian random variables.
 *
 * HISTORY:
 * Linas Vepstas October 1989
 * updates Linas Vepstas December 1992
 * new algorithms for gaussian noise, Linas Vepstas, January 1994
 */

#include <math.h>
#include "stoch.h"

#define BOX_MULLER


/* ====================================================== */

int IsStochInit = FALSE;
struct StochContext * CurrStochContext;

/* ====================================================== */

struct StochContext * CreateStochContext () {
   struct StochContext *tmp;

   tmp = (struct StochContext *) malloc (sizeof (struct StochContext));

   /* initialize default values */
#ifdef CENTRAL_LIMIT
   tmp -> niter = 30;
   tmp -> gauss_scale = sqrt (12.0 /((double) tmp -> niter));
#endif /* CENTRAL_LIMIT */

#ifdef BOX_MULLER
   tmp -> restart = TRUE;
   tmp -> next_value = 0.0;
#endif /* BOX_MULLER */

   /* initalize seed value */
   tmp -> xsubi[0] = 0;
   tmp -> xsubi[1] = 0;
   tmp -> xsubi[2] = 0;

   return (tmp);
}

/* ====================================================== */

void StochInit () {
   CurrStochContext = CreateStochContext ();
   IsStochInit = TRUE;
}

#define STOCHINIT			\
	if (!IsStochInit) StochInit ();

/* ====================================================== */

#ifdef ANSI_C
void StochSetNiter (int n)
#else
void StochSetNiter (n)
int n;
#endif
{
   STOCHINIT;

#ifdef CENTRAL_LIMIT
   CurrStochContext -> niter = n;
   /* this scale factor guarentees normal (unit std dev) distribution */
   CurrStochContext -> gauss_scale = sqrt (12.0 /((double) n));
#endif CENTRAL_LIMIT
}

/* ====================================================== */
/* 
 * This routine generates gaussian random number variables.
 * The gaussian is centered about 0.0, and has width 1.0 
 * 
 * Several algorithms are coded up here.
 *
 * #ifdef CENTRAL_LIMIT
 * The gaussian noise is generated as an (approximate) integral over
 * white noise.  The fact that gaussian noise results is due to the
 * central limit theorm.
 *
 * If X1 and X2 are two independent random variables, note the
 * varience identitiy Var (X1+X2) = Var (X1) + Var(X2)
 * Based on this identity, we can see that our normalization factor is
 * sqrt (12/N) = 2 * sqrt (3/N), the factor of 3 coming from
 * Int (-1/2 to +1/2) x**2 dx = x**3 / 3 (from -1/2 to +1/2) = 1/12
 * #endif CENTRAL_LIMIT
 *
 * The next two algorithms depend on the fact that the probability 
 * distribution of a two-dimensional gaussian is easily invertable.
 * (For a one-dimensional gaussian, we'd have to invert its distribution
 * erf(), which is hard to do. Otherwise, we could have used percentile
 * transormation to get the value directly.)
 *
 * #ifdef BOX_MULLER
 * The Box-Muller algorithm uses rejection to generate points with
 * uniform distribution within a circle, and then uses polar coordinates
 * within that circle, together percentile transormation to get gaussian
 * noise. See Athanasios Papoulis, "Probability, Random Variables, and
 * Stochasitc Processes", Third Edition, McGraw-Hill, 1991 pages 234-235.
 * #endif BOX_MULLER
 *
 * #ifdef POLAR
 * The polar-coordinate algorithm is also based on the fact that the
 * inverse of the distribution of a two-dimensional gaussian is easily
 * computable.
 * #endif POLAR
 */

double GaussianNoise () {

   STOCHINIT;

#ifdef CENTRAL_LIMIT
   {
      int i;
      double tmp;
   
      tmp = 0.0;
      for (i=0; i<CurrStochContext->niter; i++) {
         tmp += erand48 (CurrStochContext->xsubi);
      }
   
      tmp -= 0.5 * ((double) CurrStochContext->niter);
      tmp *= CurrStochContext -> gauss_scale;
      return (tmp);
   }
#endif /* CENTRAL_LIMIT */

#ifdef BOX_MULLER
   {
      double x, y, r, tmp;

      if (CurrStochContext -> restart) {
         CurrStochContext -> restart = FALSE;
         do {
            x = 2.0 * erand48 (CurrStochContext->xsubi) - 1.0;
            y = 2.0 * erand48 (CurrStochContext->xsubi) - 1.0;
            r = sqrt (x*x + y*y);
         } while (r > 1.0);
   
         tmp = sqrt (-4.0 * log(r));
         tmp /= r;
         CurrStochContext -> next_value = y * tmp;
         tmp *= x;
      } else {
         CurrStochContext -> restart = TRUE;
         tmp = CurrStochContext -> next_value;
      }
      return (tmp);
   }
#endif /* BOX_MULLER */

#ifdef POLAR
   {
      double r, phi, tmp;

      r = erand48 (CurrStochContext->xsubi);
      phi = 2.0 * M_PI * erand48 (CurrStochContext->xsubi);

      if (CurrStochContext -> restart) {
         CurrStochContext -> restart = FALSE;
         tmp = sqrt (-2.0 * log(r));
         CurrStochContext -> next_value = tmp * cos (phi);
         tmp *= sin (phi);
      } else {
         CurrStochContext -> restart = TRUE;
         tmp = CurrStochContext -> next_value;
      }
      return (tmp);
   }

#endif /* POLAR */

}

/* ====================================================== */
/*
 * This subroutine generates a white noise random variable,
 * with values between 0.0 and 1.0
 */

/*
#define WhiteNoise() drand48()
*/

/* ================= END OF FILE ========================= */
