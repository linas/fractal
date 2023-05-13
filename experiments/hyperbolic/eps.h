/*
 * eps.h 
 * Encapsulated PostScript utilities
 *
 * Linas Vepstas April 2007
 */

#ifndef __EPS_H__
#define __EPS_H__

/**
 * eps_print_prolog - set up postscript. 
 * width and height are bounding box. 
 * default size is a box from -1 to +1 with 0 at center
 */
void eps_print_prolog (int width, int height);

/**
 * Select the most basic linestyles. 
 * select line width suitable for above sizes.
 */
void eps_setup_basic_linstyles (void);

void eps_set_color_black (void);
void eps_set_color_red (void);
void eps_set_color_green (void);
void eps_set_color_blue (void);
void eps_set_color (int r, int g, int b);

/* =============================================== */

/* draw a circle of unit radius about the origin */
void eps_draw_circle(void);

/* draw a single linesegment */
void eps_draw_lineseg (double fx, double fy, double tx, double ty);

#endif
