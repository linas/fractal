/*
 * opers.h
 *
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

#ifndef OPERS_H
#define OPERS_H
/*-------------------------------------------------------------------*/

extern void absolval (float glob[],
          unsigned int sizex, unsigned int sizey);

extern void fix (float glob[],
          unsigned int sizex, unsigned int sizey);

extern double avg (float glob[],
          unsigned int sizex, unsigned int sizey);

extern double sqdev (float glob[],
          unsigned int sizex, unsigned int sizey);

extern void rescale (float glob[],
             unsigned int sizex, unsigned int sizey,
             float scale_factor);

extern void expand (float glob[],
             unsigned int sizex, unsigned int sizey,
             float scale_factor, float offset);

extern void recip (float glob[],
              unsigned int sizex, unsigned int sizey);

extern void takelog (float glob[],
              unsigned int sizex, unsigned int sizey);

extern void fakelog (float glob[],
              unsigned int sizex, unsigned int sizey,
              float scale_factor);

extern void clamp (float glob[],
             unsigned int sizex, unsigned int sizey,
             float scale_factor, float offset);

extern void dump (float glob[],
             unsigned int sizex, unsigned int sizey);

#endif  /* OPERS_H */

/* --------------- END OF FILE ------------------------- */
