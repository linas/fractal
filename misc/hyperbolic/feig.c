/*
 * feig.c
 *
 * Explore possible feignbaum ratios in the forking of the tree.
 * Derived as a heavily-edited version of tree.c
 * Conclusion -- the feigenbaum ratio here seems to converge
 * to the limit of 1.0 ... (!?)  The feigenbaum const does NOT
 * appear naturally.  Not too surprising ... things would be too
 * easy if it did.  It does beg the question, though: there is, after
 * all, something fundamentall different between the ordinary
 * period-doubling maps, and the hyperbolic plane. What, then, is the
 * difference?  What's really going on with the maps?
 *
 * Linas Vepstas August 2016
 */

#include "cplex.h"
#include "flt.h"
#include "question.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double pnm2 = 1.0;
double pnm1 = 0.5;
double pn = 0.25;

double qnm2 = 1.0;
double qnm1 = 0.5;
double qn = 0.25;

/* arc */
void print_arc(mobius_t m, cplex zstart, cplex zend)
{
	cplex za = mobius_xform (m, zstart);
	cplex zb = mobius_xform (m, zend);
	double dist = cplex_dist (za, zb);
	printf("duuude re=%Lg im=%Lg dist=%g\n", za.re, za.im, dist);
	// printf("other end= %Lg %Lg\n", zb.re, zb.im);
	printf("basic ratio= %Lg %Lg\n", za.re/zb.re, za.im/zb.im);

	pnm2 = pnm1;
	pnm1 = pn;
	pn = za.im;
/*
	pn = dist;
	pn = za.im * za.im;
	pn = dist * dist;
	pn = sqrt(za.im);
	pn = sqrt(dist);
	pn = log(za.im);
	// pn = question_mark (i,deno);
	pn = question_inverse (za.im);
*/

	double feig = (pnm1-pnm2)/(pn-pnm1);
	printf("feig= %g %g\n", feig, 4.0*feig);

	qnm2 = qnm1;
	qnm1 = qn;
	qn = feig;
	double rat = (qnm1-qnm2)/(qn-qnm1);

	// printf("frat= %g %g\n", rat, 4.0*rat);
}

/* ==================================================== */

mobius_t go_to_fork_tip(double sign)
{
	cplex c = cplex_set (-0.5, 0.0);
	mobius_t m = disk_center (c);
	cplex rot = cplex_exp_itheta (sign*2.0*M_PI/6.0);
	m = mobius_scale (m,rot);
// show_mobius (m);
	return m;
}

void draw_fork(mobius_t m, int level)
{
	if (level == 0) return;
	level--;
printf("----- level=%d\n", level);

	cplex za, zb;

	za = cplex_set (0.0, 0.0);

	zb = cplex_set (0.25, 0.25*sqrt(3.0));
	// za = mobius_xform (m,za);
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	print_arc (m, za, zb);

/****
	zb = cplex_set (0.25, -0.25*sqrt(3.0));
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	print_arc (m, za, zb);
****/

	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);

/***
	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);
***/
}


void draw()
{
	mobius_t m;
	int level=3;
level=150;

	double cent = sqrt(3.0) - 2.0;
	cplex z = cplex_set (cent, 0.0);
	// z = cplex_set (-0.25, 0.0); xxxxxx
	z = cplex_set (0.0, 0.0);
	mobius_t off = disk_center (z);

	mobius_t rot = mobius_rotate (-0.5*M_PI);
	// mobius_t rot = mobius_rotate ((-0.5-0.166666)*M_PI);
	off = mobius_mul (rot, off);

#define HALF_PLANE
#ifdef HALF_PLANE
	mobius_t xfm = mobius_set (1,0,0,1);
	xfm = to_half_plane_xform();
	off = mobius_mul (xfm, off);
#endif

	m = go_to_fork_tip(-1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{
	draw();
}
