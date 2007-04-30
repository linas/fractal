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

	cplex tip = cplex_set(0.5,0.5*sqrt(3.0));

eps_set_color_red();
	mobius_t tee = mobius_set (1,1,0,1);
	cplex zb = mobius_xform (tee, tip);
	draw_seg (m, tip, zb);

	mobius_t mt = mobius_mul (m,tee);
	recursive_draw (depth, mt);

eps_set_color_green();
	mobius_t ess = mobius_set (0,-1,1,0);
	zb = mobius_xform (ess, tip);
	draw_seg (m, tip, zb);

	mobius_t ms = mobius_mul (m,ess);
	recursive_draw (depth, ms);
}

void draw (int n)
{
	mobius_t ident = mobius_ident();

	recursive_draw (16,ident);
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
