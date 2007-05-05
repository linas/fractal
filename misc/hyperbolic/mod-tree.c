/*
 * tree.c 
 * draw a binary tree on the hyperbolic plane
 * cosntructed from modular xforms
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include "eps.h"
#include "flt.h"
#include "flt-eps.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void eps_setup_plane (void)
{
   eps_setup_basic_linstyles();
	printf ("0.0 -0.8 translate\n");
	printf ("0.4 1.6 scale\n");
}

void recursive_draw (int depth, mobius_t m)
{
	if (0 >= depth)
		return;
	depth --;

	// cplex tip = cplex_set(0.5,0.5*sqrt(3.0));
	cplex tip = cplex_set(-0.5,0.5*sqrt(3.0));

eps_set_color_red();
	// mobius_t are = mobius_set (1,1,0,1);
	mobius_t are = mobius_set (1,-1,0,1);
	cplex zb = mobius_xform (are, tip);
	draw_seg (m, tip, zb);

	mobius_t mr = mobius_mul (m,are);
	recursive_draw (depth, mr);

eps_set_color_green();
	// mobius_t ell = mobius_set (1,0,1,1);
	mobius_t ell = mobius_set (1,0,-1,1);
	zb = mobius_xform (ell, tip);
	draw_seg (m, tip, zb);

	mobius_t ml = mobius_mul (m,ell);
	recursive_draw (depth, ml);
}

void draw (int n)
{
	mobius_t ident = mobius_ident();

	recursive_draw (n,ident);
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{
	if (argc < 1)
	{
		fprintf (stderr, "Usage: %s <shift>\n", argv[0]);
		exit(1);
	}

	int n = atoi (argv[1]);

	eps_print_prolog(800,200);
	eps_setup_plane();

	draw(n);

	printf ("showpage\n");
}
