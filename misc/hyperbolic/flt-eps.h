/*
 * flt-eps.h 
 * fractional linear transfor (flt) postscript drawing utils
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include "flt.h"

/* Draw straight line segment */
void draw_seg(mobius_t m, cplex zf, cplex zt);

/* draw line segment in the klein model */
void draw_klein_seg(mobius_t m, cplex zf, cplex zt);

/* draw statically-tesselated arc */
void draw_tesselated_arc(mobius_t m, cplex zf, cplex zt);

/* draw dynamically-tesselated arc */
void draw_arc(mobius_t m, cplex zf, cplex zt);

