/*
 * flt.h 
 * Generic Fractional Linear Transform (Mobius transform) 
 * definition and handling routines
 *
 * Linas Vepstas April 2007
 */

#ifndef __LFUNC_FLT_H__
#define __LFUNC_FLT_H__

#include <math.h>
#include <stdio.h>
#include "cplex.h"

/* fractional linear transform */
typedef struct {
	cplex a,b,c,d;
} mobius_t;

/* Create an element of SL(2,R). Its up to user to ensure
 * that ad-bc=1 */
static inline mobius_t mobius_set(long double a, long double b, 
                                  long double c, long double d)
{
	mobius_t m;
	m.a = cplex_set (a, 0.0L);
	m.b = cplex_set (b, 0.0L);
	m.c = cplex_set (c, 0.0L);
	m.d = cplex_set (d, 0.0L);
	return m;
}

/* Return the identity matrix */
static inline mobius_t mobius_ident(void)
{
	return mobius_set (1.0L,0.0L,0.0L,1.0L);
}

/* Compute teh determinant */
static inline cplex mobius_determinant(mobius_t m)
{
	return cplex_sub (cplex_mul (m.a, m.d), cplex_mul(m.b, m.c));
}

/* Return the inverse */
static inline mobius_t mobius_invert(mobius_t m)
{
	cplex det = mobius_determinant(m);
	det = cplex_recip(det);

	mobius_t n;
	n.a = cplex_mul(m.d, det);
	n.d = cplex_mul(m.a, det);
	n.b = cplex_mul(cplex_neg(m.b), det);
	n.c = cplex_mul(cplex_neg(m.c), det);
	return n;
}

/* Return a transformation that simply rotates by theta radians */
static inline mobius_t mobius_rotate (long double theta)
{
	mobius_t m;
	m.a = cplex_exp_itheta (theta);
	m.b = cplex_zero();
	m.c = cplex_zero();
	m.d = cplex_one();
	return m;
}

/* Perform a multiplicative scaling of the thing */
static inline mobius_t mobius_scale(mobius_t m, const cplex z)
{
	m.a = cplex_mul (m.a, z);
	m.b = cplex_mul (m.b, z);
	return m;
}

/* Product of mobius transforms (just a matrix multiply) */
static inline mobius_t mobius_mul(const mobius_t l, const mobius_t r)
{
	mobius_t m;

	m.a = cplex_add (cplex_mul (l.a, r.a), cplex_mul(l.b, r.c));
	m.c = cplex_add (cplex_mul (l.c, r.a), cplex_mul(l.d, r.c));
	m.b = cplex_add (cplex_mul (l.a, r.b), cplex_mul(l.b, r.d));
	m.d = cplex_add (cplex_mul (l.c, r.b), cplex_mul(l.d, r.d));

	return m;
}

/* apply mobius xform to z, return (az+b)/(cz+d) */
static inline cplex mobius_xform (const mobius_t m, const cplex z)
{
	cplex numer = cplex_mul(m.a, z);
	numer = cplex_add (numer, m.b);
	cplex deno = cplex_mul(m.c, z);
	deno = cplex_add (deno, m.d);
	cplex w = cplex_div (numer,deno);
	return w;
}

/* Return mobius xform that recenters the disk at w */
static inline mobius_t disk_center (cplex w)
{
	cplex nu = w;
	nu.re += 1.0L;
	cplex de = w;
	de.re -= 1.0L;

	cplex z = cplex_div (nu,de);
	z = cplex_neg (z);
	z = cplex_times_i (z);
	cplex zb = cplex_conj (z);

	cplex mi = cplex_set(0.0L,-1.0L);

	mobius_t m;
	m.a = cplex_sub (mi, z);
	m.b = cplex_add (mi, z);
	m.c = cplex_sub (mi, zb);
	m.d = cplex_add (mi, zb);
	return m;
}

/**
 * Return a mobius transform that maps the point z=0 (the 
 * center of the Poincare disk) to z=i of the half-plane.
 * More generally, the frac lin xform will take any
 * point of the disk and map it to a point in the half-plane.
 */
static inline mobius_t to_half_plane_xform(void)
{
	mobius_t m;
	m.a = cplex_set (0.0L, -1.0L);
	m.b = cplex_set (0.0L, -1.0L);
	m.c = cplex_set (1.0L, 0.0L);
	m.d = cplex_set (-1.0L, 0.0L);
	return m;
}

/**
 * Return a mobius transform that maps the point z=+i to the 
 * center of the Poincare disk */
static inline mobius_t to_disk_xform(void)
{
	mobius_t m;
#if XX
	m.a = cplex_set (1.0L, 0.0L);
	m.b = cplex_set (0.0L, -1.0L);
	m.c = cplex_set (1.0L, 0.0L);
	m.d = cplex_set (0.0L, 1.0L);
#endif

	m.a = cplex_set (0.0L, 0.5L);
	m.b = cplex_set (0.5L, 0.0L);
	m.c = cplex_set (0.0L, 0.5L);
	m.d = cplex_set (-0.5L, 0.0L);
	return m;
}

/**
 * Convert a half-plane xform to a disk xform 
 * This is just a similarity xform given by
 * return mobius_mul( mobius_mul (to_disk_xform(), m), to_half_plane_xform())
 */
static inline mobius_t to_disk(mobius_t m)
{
	cplex apd = cplex_add (m.a, m.d);
	cplex amd = cplex_sub (m.a, m.d);
	cplex bpc = cplex_times_i (cplex_add (m.b, m.c));
	cplex bmc = cplex_times_i (cplex_sub (m.b, m.c));

	mobius_t n;
	n.a = cplex_scale (0.5L, cplex_add (apd, bmc));
	n.b = cplex_scale (0.5L, cplex_sub (amd, bpc));
	n.c = cplex_scale (0.5L, cplex_add (amd, bpc));
	n.d = cplex_scale (0.5L, cplex_sub (apd, bmc));
	
	return n;
}

/**
 * Convert a disk xform to a plane xform 
 * This is just a similarity xform given by
 * return mobius_mul( mobius_mul (to_half_plane_xform(), m), to_disk_xform())
 */
static inline mobius_t to_half_plane(mobius_t m)
{
	cplex apd = cplex_add (m.a, m.d);
	cplex amd = cplex_sub (m.a, m.d);
	cplex bpc = cplex_add (m.b, m.c);
	cplex bmc = cplex_sub (m.b, m.c);

	mobius_t n;
	n.a = cplex_scale (0.5L, cplex_add (apd, bpc));
	n.b = cplex_scale (0.5L, cplex_times_i (cplex_sub (amd, bmc)));
	n.b = cplex_neg (n.b);
	n.c = cplex_scale (0.5L, cplex_times_i (cplex_add (amd, bmc)));
	n.d = cplex_scale (0.5L, cplex_sub (apd, bpc));
	
	return n;
}

static inline void show_mobius(mobius_t m)
{
	printf ("a=%Lg %+Lgi    b=%Lg %+Lgi\n", m.a.re, m.a.im, m.b.re, m.b.im);
	printf ("c=%Lg %+Lgi    d=%Lg %+Lgi\n", m.c.re, m.c.im, m.d.re, m.d.im);
}

#endif /* __LFUNC_FLT_H__ */

/* ======================= END OF FILE ============================= */

