/* 
 * NAME:
 * image.h
 *
 * FUNCTION:
 * generate pixmaps showing Poincare recurrence time for circle map.
 *
 * HISTORY:
 * created Linas Vepstas Oct 1989
 * major revs Linas Vepstas Jan 1991
 * spruced up Linas Vepstas July 1993
 * fixed major performance bug (AIX fmod() stinks!) -- January 1994
 * added logistic map -- February 1994
 * added classic mandelbrot -- June 1995
 */

typedef float (*DensityCB) (void*, double, double);

/*-------------------------------------------------------------------*/
/**
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 */

void walk_rect (float *glob, 
                unsigned int sizex,  /* width, in pixels */
                unsigned int sizey,  /* height, in pixels */
                double x_min,        /* left side of pixmap */
                double x_max,        /* right side of pixmap */
                double y_min,        /* bottom of pixmap */
                double y_max,        /* top of pixmap */
                float (*callback)(void *, double, double),  /* callback */
                void * calldata);     /* static data */
/*-------------------------------------------------------------------*/
/** 
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 * instead of filling in in a linear fashion, its done with a
 * trianglelar distortion.
 */

void walk_tri (float *glob, 
               unsigned int sizex,  /* width, in pixels */
               unsigned int sizey,  /* height, in pixels */
               double x_min,        /* left side of pixmap */
               double x_max,        /* right side of pixmap */
               double y_min,        /* bottom of pixmap */
               double y_max,        /* top of pixmap */
               float (*callback)(void *, double, double),  /* callback */
               void * calldata);     /* static data */
   
/*-------------------------------------------------------------------*/
/* 
 * This routine visits all the pixels of a pixmap, filling them in with
 * values returned by the supplied callback.
 * instead of filling in in a linear fashion, its done with a
 * trianglelar distortion.
 */

void walk_utri (float *glob, 
                unsigned int sizex,  /* width, in pixels */
                unsigned int sizey,  /* height, in pixels */
                double x_min,        /* left side of pixmap */
                double x_max,        /* right side of pixmap */
                double y_min,        /* bottom of pixmap */
                double y_max,        /* top of pixmap */
                float (*callback)(void *, double, double),  /* callback */
                void * calldata);     /* static data */

/*-------------------------------------------------------------------*/
