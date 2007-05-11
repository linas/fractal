/*
 * flt-eps.c 
 * fractional linear transfor (flt) postscript drawing utils
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include "eps.h"
#include "flt.h"
#include <math.h>
#include <stdio.h>

/* ==================================================== */

/* Draw straight line segment */
void draw_seg(mobius_t m, cplex zf, cplex zt)
{
	cplex za = mobius_xform (m, zf);
	cplex zb = mobius_xform (m, zt);

	// clip anything taller than yclip
	double yclip = 5;
	if (yclip <za.im)
	{
		double m = (zb.im-za.im)/(zb.re-za.re);
		double xclip = zb.re +(yclip-zb.im)/m;
		za.re = xclip;
		za.im = yclip;
	}
	if (yclip <zb.im)
	{
		double m = (zb.im-za.im)/(zb.re-za.re);
		double xclip = za.re +(yclip-za.im)/m;
		zb.re = xclip;
		zb.im = yclip;
	}
	eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
}

/* Draw line segment in the klein model */
void draw_klein_seg(mobius_t m, cplex zf, cplex zt)
{
	cplex za = mobius_xform (m, zf);
	cplex zb = mobius_xform (m, zt);

	double ma = 2.0/(1.0 + za.re * za.re + za.im * za.im);
	double mb = 2.0/(1.0 + zb.re * zb.re + zb.im * zb.im);
	za = cplex_scale (ma,za);
	zb = cplex_scale (mb,zb);
	eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
}

/* Draw statically-tesselated arc */
void draw_tesselated_arc(mobius_t m, cplex zf, cplex zt)
{
	int i;
	// int nseg=20;
	int nseg=5;

	cplex zstart = zf;
	cplex zdelta = cplex_scale((1.0/(double)nseg), cplex_sub (zt,zf));
	cplex za = mobius_xform (m, zstart);

	for (i=0; i<nseg; i++)
	{
		cplex zend = cplex_add (zstart, zdelta);
		cplex zb = mobius_xform (m, zend);
		eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
		zstart = zend;
		za = zb;
	}
}

#if 0
/* do the klein geometry, where all arcs are straight lines */
void draw_arc(mobius_t m, cplex zf, cplex zt)
{
	draw_klein_seg(m,zf,zt);
}
#endif

/* draw dynamically-tesselated arc */
void draw_arc(mobius_t m, cplex zf, cplex zt)
{
	int i;
	int nseg=20;

	// clip anything taller than yclip
	double yclip = 5;
	if (yclip <zf.im)
	{
		double m = (zt.im-zf.im)/(zt.re-zf.re);
		double xclip = zt.re +(yclip-zt.im)/m;
		zf.re = xclip;
		zf.im = yclip;
	}
	if (yclip <zt.im)
	{
		double m = (zt.im-zf.im)/(zt.re-zf.re);
		double xclip = zf.re +(yclip-zf.im)/m;
		zt.re = xclip;
		zt.im = yclip;
	}

	cplex zstart = zf;
	cplex zdelta = cplex_scale((1.0/(double)nseg), cplex_sub (zt,zf));
	cplex za = mobius_xform (m, zstart);
	cplex zb;

	cplex zs = za;

	for (i=0; i<nseg; i++)
	{
		cplex zend = cplex_add (zstart, zdelta);
		zb = mobius_xform (m, zend);
		double dist = cplex_dist (zs, zb);
		if (dist > 0.03) {
			eps_draw_lineseg (zs.re, zs.im, zb.re, zb.im);
			zs = zb;
		}
		zstart = zend;
		za = zb;
	}

	double dist = cplex_dist (zs, zb);
	if (dist > 0.001) 
		eps_draw_lineseg (zs.re, zs.im, zb.re, zb.im);
}

/* ================== END OF FILE ============= */
