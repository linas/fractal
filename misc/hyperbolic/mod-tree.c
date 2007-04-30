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

void eps_setup_plane (void)
{
   eps_setup_basic_linstyles();
	printf ("0.0 -0.8 translate\n");
	printf ("0.4 1.6 scale\n");
}

void draw (int n)
{
	mobius_t ident = mobius_ident();
	mobius_t tee = mobius_set (1,1,0,1);

	cplex eye = cplex_set(0,1);

	cplex zb = mobius_xform (tee, eye);

	draw_seg (ident, eye, zb);
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{
	int n = atoi (argv[1]);

	eps_print_prolog(400,100);
	eps_setup_plane();

	draw(n);

	printf ("showpage\n");
}
