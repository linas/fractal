
#ifdef AIX221
#ifndef _C_func
#define _C_func
#endif
#endif

#ifdef AIX315
#define ANSI_C
#endif

#include <math.h>


#ifdef ANSI_C
#include <stddefs.h>
#else
#include <sys/types.h>
#endif /* ANSI_C */

#ifdef ANSI_C
extern double drand48 (void);
extern double erand48 (unsigned short xsubi[3]);
#else
extern double drand48 ();
extern double erand48 ();
#endif /* ANSI_C */

extern int IsStochInit;

struct StochContext {

   /* the following two are private and are used only by the 
    * central limit algorithm */
   int		niter;
   double	gauss_scale;

   /* the following are private and are used only by the 
    * box-muller algorithm */
   int		restart;
   double	next_value;

   /* random number generator context for erand48() */
   unsigned short xsubi[3];
};


/*
 * won't do this, since we do not know if we initalized 
 * #define WhiteNoise() erand48(CurrStochContext->xsubi)
 */

#define WhiteNoise() drand48()

/* ================================================= */

#ifdef ANSI_C

extern void StochInit (void);
extern void StochSetNiter (int n); /* this routine used ONLY by 
                                    * central-limit algorithm */
extern double GaussianNoise (void);

#else

extern void StochInit ();
extern void StochSetNiter ();
extern double GaussianNoise ();

#endif

/* ================= END OF FILE ====================== */
