/*
 * opers.c
 * 
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */
#include <stdio.h>
#include <math.h>

/*-------------------------------------------------------------------*/

void fix (float glob[], 
          unsigned int sizex, unsigned int sizey)
{
   long		i;
   for (i=0; i<sizex*sizey; i++) {
      if (glob[i] <= 0.0) glob [i] = 1.0/20000.0;
   }
}

/*-------------------------------------------------------------------*/

void absolval (float glob[], 
          unsigned int sizex, unsigned int sizey)
{
   long		i;
   for (i=0; i<sizex*sizey; i++) {
      glob[i] = fabs (glob[i]);
   }
}

/*-------------------------------------------------------------------*/

double avg (float glob [], 
            unsigned int sizex, unsigned int sizey)
{
   long		i;
   double retval;

   /* renormalize */
   retval = 0.0;
   for (i=0; i<sizex*sizey; i++) {
      retval += glob [i];
   }
   retval /= (double) sizex*sizey;
   return retval;
}

/*-------------------------------------------------------------------*/

double sqdev (float glob [], 
              unsigned int sizex, unsigned int sizey)
{
   long		i;
   double avg, msq;

   /* renormalize */
   avg = 0.0;
   msq = 0.0;
   for (i=0; i<sizex*sizey; i++) {
      avg += glob [i];
      msq += (glob[i]*glob [i]);
   }
   
   avg /= (double) sizex*sizey;
   msq /= (double) sizex*sizey;
   msq -= avg*avg;
   return msq;
}

/*-------------------------------------------------------------------*/

void rescale (float glob[], 
              unsigned int sizex, unsigned int sizey, 
              float scale_factor)
{
   long		i;
   int n;

   /* renormalize */
   for (i=0; i<sizex*sizey; i++) {
      glob [i] *= scale_factor ;
      n = (int) glob[i];
      if (0.0 > glob[i]) n--;
      glob [i] -= n;
   }
}

/*-------------------------------------------------------------------*/

void expand (float glob[], 
             unsigned int sizex, unsigned int sizey, 
             float scale_factor, float offset)
{
   long		i;
   
   /* renormalize */
   for (i=0; i<sizex*sizey; i++) {
      glob [i] += offset;
      glob [i] *= scale_factor ;
   }
}

/*-------------------------------------------------------------------*/

void takelog (float glob[], 
              unsigned int sizex, unsigned int sizey)
{
   long		i;
   
   /* renormalize */
   for (i=0; i<sizex*sizey; i++) {
      glob [i] = (float) log ((double) glob[i]);
   }
}

/*-------------------------------------------------------------------*/

void fakelog (float glob[], 
              unsigned int sizex, unsigned int sizey, 
              float scale_factor)
{
   long		i;
   float		ska;
   
   ska = exp (scale_factor) - 1.0;
   /* renormalize */
   for (i=0; i<sizex*sizey; i++) {
      glob [i] = (float) log ((double) (ska * glob[i] + 1.0)) / scale_factor;
   }
}

/*-------------------------------------------------------------------*/

void clamp (float glob[], 
              unsigned int sizex, unsigned int sizey, 
              float clmin, float clmax)
{
   long		i;

   /* clamp */
   for (i=0; i<sizex*sizey; i++) {
      if (glob[i] < clmin) glob[i] = clmin;
      if (glob[i] > clmax) glob[i] = clmax;
   }
}

/*-------------------------------------------------------------------*/

void dump (float glob[], 
              unsigned int sizex, unsigned int sizey)
{
   long		i;
   
   /* renormalize */
   for (i=0; i<sizex*sizey; i++) {
      printf ("(%ld, %ld) = %g \n", i%sizex, i/sizex,  glob[i]);
   }
}

/* --------------- END OF FILE ------------------------- */
