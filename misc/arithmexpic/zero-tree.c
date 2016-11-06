/*
 * Find zeros of function.
 *
 * Assumes that the zeros can be organized into a binary tree,
 * and that the tree is the Stern-Brocot tree, so that the zeros
 * are located according to the tree.
 *
 * April 2016, October 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>
#include <mp-complex.h>
#include <mp-consts.h>
#include <mp-misc.h>
#include <mp-trig.h>
#include <mp-zerofind.h>

#include "moebius.h"

#include "genfunc.h"

// Partition function
void parti_z_mpf(mpf_t res, long n)
{
	mpz_t part; mpz_init(part);
	partition_z(part, n);
	mpf_set_z(res, part);
	mpz_clear(part);
}

void divisor_mpf(mpf_t res, long n)
{
	mpf_set_ui(res, divisor(n));
}

void parti(cpx_t f, cpx_t z, int nprec)
{
	// cpx_exponential_genfunc_mpf(f, z, nprec, parti_z_mpf);
	cpx_exponential_genfunc_mpf(f, z, nprec, divisor_mpf);
}

bool survey_cell(void (*func)(cpx_t f, cpx_t z, int nprec),
                 double rguess, double tguess, double cell_size,
                 int ndigits, int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;

	cpx_t a,b,c,d;
	cpx_init2(a, bits);
	cpx_init2(b, bits);
	cpx_init2(c, bits);
	cpx_init2(d, bits);

	cpx_t fa,fb,fc,fd;
	cpx_init2(fa, bits);
	cpx_init2(fb, bits);
	cpx_init2(fc, bits);
	cpx_init2(fd, bits);

	// Set up four corners
	double st = sin(tguess);
	double ct = cos(tguess);
	cpx_set_d(a, rguess*ct, rguess*st);
	double rgd = rguess + cell_size;
	cpx_set_d(b, rgd*ct, rgd*st);
	tguess += cell_size / rguess;
	st = sin(tguess);
	ct = cos(tguess);
	cpx_set_d(c, rgd*ct, rgd*st);
	cpx_set_d(d, rguess*ct, rguess*st);

	// Evaluate at the four corners
	func(fa, a, nprec);
	func(fb, b, nprec);
	func(fc, c, nprec);
	func(fd, d, nprec);

	// Compute contour integral i.e. sum the phases.
	double phase = atan2(cpx_get_im(fa), cpx_get_re(fa));
	double sum = phase;
	double hi = phase;
	double lo = phase;

	phase = atan2(cpx_get_im(fb), cpx_get_re(fb));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	phase = atan2(cpx_get_im(fc), cpx_get_re(fc));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	phase = atan2(cpx_get_im(fd), cpx_get_re(fd));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	// We expect the phase to wind around, so that hi and low
	// differ by almost 2pi.
	// We expect the integral to be large, approaching 2pi.
	// if (hi-lo < 4.0 or fabs(phase) < 3.0)
	if (hi-lo < 4.0)
	{
		cpx_clear(a);
		cpx_clear(b);
		cpx_clear(c);
		cpx_clear(d);

		cpx_clear(fa);
		cpx_clear(fb);
		cpx_clear(fc);
		cpx_clear(fd);

		return false;
	}

	// Now, a will be the center of the rectangle.
	cpx_add(a, a, b);
	cpx_add(a, a, c);
	cpx_add(a, a, d);
	cpx_div_ui(a, a, 4);

	cpx_t e1, e2, zero;
	cpx_init2(e1, bits);
	cpx_init2(e2, bits);
	cpx_init2(zero, bits);
	cpx_sub(e1, b, a);
	cpx_sub(e2, d, a);

	int rc = cpx_find_zero(zero, func, a, e1, e2, ndigits, nprec);

if (rc) {printf("duuude found noothing\n");}

	// if rc is not zero, then nothing was found
	if (rc)
	{
		cpx_clear(a);
		cpx_clear(b);
		cpx_clear(c);
		cpx_clear(d);

		cpx_clear(fa);
		cpx_clear(fb);
		cpx_clear(fc);
		cpx_clear(fd);

		cpx_clear(e1);
		cpx_clear(e2);
		return false;
	}

	printf("---------------\n");
	cpx_prt("zero = ", zero); printf("\n");

	mpf_t r, t, pi;
	mpf_init2(r, bits);
	mpf_init2(t, bits);
	mpf_init2(pi, bits);

	cpx_abs(r, zero);
	fp_prt("r = ", r); printf("\n");

	fp_arctan2(t, zero->im, zero->re, nprec);
	fp_pi(pi, nprec);
	mpf_div(t, t, pi);

	fp_prt("theta/pi = ", t); printf("\n");

#define CHECK_RESULT 1
#ifdef CHECK_RESULT
	cpx_t check;
	cpx_init(check);
	func(check, zero, nprec);
	double eps_r = cpx_get_re(check);
	double eps_i = cpx_get_im(check);
	double eps = sqrt(eps_r*eps_r + eps_i*eps_i);
	printf("fun at zero = %g\n", eps);
#endif

	printf("\n");
	fflush(stdout);

	cpx_clear(a);
	cpx_clear(b);
	cpx_clear(c);
	cpx_clear(d);

	cpx_clear(fa);
	cpx_clear(fb);
	cpx_clear(fc);
	cpx_clear(fd);

	cpx_clear(e1);
	cpx_clear(e2);
	return true;
}

void survey(void (*func)(cpx_t f, cpx_t z, int nprec),
            double rmax, double cell_size,
            int ndigits, int nprec)
{
	for (double r=1.0; r<rmax; r += cell_size)
	{
		double step = cell_size / r;
		for (double t=0.0; t < M_PI; t += step)
		{
			survey_cell(func, r, t, cell_size, ndigits, nprec);
		}
	}
}

// =================================================================
// Ordinary fraction numer/denom
typedef struct fraction
{
	long numer;
	long denom;
} Fraction;

// Tree location, expressed several ways: as dyadic, Farey, and bounds.
// This is only a symbolic labelling of the tree node.
typedef struct tree_node
{
	// Index runs from 0 to infinity
	int index;

	// index = 2^level + vindex
	// level runs from 0 to infty, level = number of bits in index.
	int level;
	// vindex runs from 0 to 2^(level-1)
	int vindex;

	// denominator of dyadic fraction is 2^level
	// numerator of dyadic fraction is 2*vindex +1
	Fraction dyadic;

	// Farey fraction corresponding to above.
	Fraction farey;

	// Left and right parents of the Farey fraction
	Fraction left;
	Fraction right;
} TreeNode;

// Initialize the node to the indicated location.
void make_node(int idx, TreeNode* node)
{
	node->index = idx;

	// get the level and offset in the level
	int lvl = 0;
	while (idx) { lvl++; idx >>= 1;}
	node->level = lvl;

	idx = node->index;
	int dyden = 1 << lvl;
	node->vindex = idx - dyden/2;

	// Set up the dyadic fraction; this is straight-forward
	node->dyadic.numer = 2 * node->vindex +1;
	node->dyadic.denom = dyden;

	// Get the Farey fractions recursively.
	if (0 == idx)
	{
		node->left.numer = 0;
		node->left.denom = 1;
		node->right.numer = 1;
		node->right.denom = 1;
		node->farey.numer = 0;
		node->farey.denom = 1;
		return;
	}

	TreeNode rec;
	make_node(idx/2, &rec);
	if (0 == idx%2)
	{
		node->left = rec.left;
		node->right = rec.farey;
	}
	else
	{
		node->left = rec.farey;
		node->right = rec.right;
	}

	node->farey.numer = node->left.numer + node->right.numer;
	node->farey.denom = node->left.denom + node->right.denom;
}

// ==================================================================

// Numeric description of the zero, and the expeted bounds for the zero.
typedef struct zero_node
{
	TreeNode tree_node;
	double phi_min;
	double phi_max;
	mpf_t  phi_farey;
} ZeroNode;

void make_zero(int idx, ZeroNode* zn)
{
	make_node(idx, &zn->tree_node);

	// Set up the left and right bounds.
	long n = zn->tree_node.left.numer;
	long d = zn->tree_node.left.denom;
	zn->phi_min = 2.0 * ((double) n) / ((double) d);
	n = zn->tree_node.right.numer;
	d = zn->tree_node.right.denom;
	zn->phi_max = 2.0 * ((double) n) / ((double) d);

	// Set up the expected location of the zero.
	mpf_init(zn->phi_farey);
	mpf_set_ui(zn->phi_farey, 2 * zn->tree_node.farey.numer);
	mpf_div_ui(zn->phi_farey, zn->phi_farey, zn->tree_node.farey.denom);
}

void clear_zero(ZeroNode* zn)
{
	mpf_clear(zn->phi_farey);
}

// ==================================================================

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <tree-depth>\n", argv[0]);
		exit(1);
	}

	int nprec = 80;
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;
	mpf_set_default_prec(bits);

	int lvl = atoi(argv[1]);
	int midx = 1<<lvl;

	for (int idx=0; idx < midx; idx++)
	{
		// Skip the right side of the tree.
		int tb = idx;
		while (3 <= tb) { tb >>= 1; }
		if (2 > tb) continue;

		ZeroNode zero;
		make_zero(idx, &zero);

		TreeNode node;
		// make_node(idx, &node);
		node = zero.tree_node;

		printf("node ");
		// printf("idx=%d, lv=(%d, %d) ", idx, node.level, node.vindex);
		printf("dy=%ld/%ld ", node.dyadic.numer, node.dyadic.denom);
		printf("min=%ld/%ld ", node.left.numer, node.left.denom);
		printf("max=%ld/%ld ", node.right.numer, node.right.denom);
		printf("far=%ld/%ld ", node.farey.numer, node.farey.denom);
		printf("phi = (%g %g) ", zero.phi_min, zero.phi_max);
		printf("phx= %g ", mpf_get_d(zero.phi_farey));
		printf("\n");

		clear_zero(&zero);
	}

	printf("Done!\n");
}
