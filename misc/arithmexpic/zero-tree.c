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
                 cpx_t zero,
                 double rguess, double tguess,
                 double rdelta, double tdelta,
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

	tguess *= M_PI;
	tdelta *= M_PI;

	// Set up four corners, on each side of the guess center.
	// a,b,c,d come in right-handed order.
	double st = sin(tguess-tdelta);
	double ct = cos(tguess-tdelta);
	double rgd = rguess - rdelta;
	cpx_set_d(a, rgd*ct, rgd*st);

	rgd = rguess + rdelta;
	cpx_set_d(b, rgd*ct, rgd*st);

	st = sin(tguess+tdelta);
	ct = cos(tguess+tdelta);
	cpx_set_d(c, rgd*ct, rgd*st);

	rgd = rguess - rdelta;
	cpx_set_d(d, rgd*ct, rgd*st);

	// Evaluate at the four corners
	func(fa, a, nprec);
	func(fb, b, nprec);
	func(fc, c, nprec);
	func(fd, d, nprec);

	// Compute contour integral i.e. sum the phases.
	double phase = atan2(cpx_get_im(fa), cpx_get_re(fa));
	if (phase < 0.0) phase += 2.0 * M_PI;
	double sum = phase;
	double hi = phase;
	double lo = phase;

	phase = atan2(cpx_get_im(fb), cpx_get_re(fb));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

	phase = atan2(cpx_get_im(fc), cpx_get_re(fc));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

	phase = atan2(cpx_get_im(fd), cpx_get_re(fd));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

// printf("duuude explore at r=%g t=%g delta=%g sum=%g\n", rguess, tguess/M_PI, hi-lo, 0.25*sum);
	// We expect the phase to wind around, so that hi and low
	// differ by almost 2pi.
	// We expect the integral to be large, approaching pi,
	// but if we are unlucky, as low as 3pi/8, I guess.
	if (hi-lo < 4.0 && 0.25*sum < 1.5)
	{
		cpx_clear(a);
		cpx_clear(b);
		cpx_clear(c);
		cpx_clear(d);

		cpx_clear(fa);
		cpx_clear(fb);
		cpx_clear(fc);
		cpx_clear(fd);

// printf("duuude reject at r=%g t=%g hi=%g lo=%g sum=%g\n", rguess, tguess/M_PI, hi, lo, 0.25*sum);
		return false;
	}

// printf("duude candidate at %g %g\n", rguess, tguess/M_PI);
	// Now, set a to be the best-guess; its in the center of the rectangle.
	st = sin(tguess);
	ct = cos(tguess);
	cpx_set_d(a, rguess*ct, rguess*st);

	cpx_t e1, e2;
	cpx_init2(e1, bits);
	cpx_init2(e2, bits);
	cpx_sub(e1, b, a);
	cpx_sub(e2, d, a);

	int rc = cpx_find_zero(zero, func, a, e1, e2, ndigits, nprec);

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
// printf("duuude found nothing\n");
		return false;
	}

	// cpx_prt("zero = ", zero); printf("\n");

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

	// We expect the zero to be found between these two bounds.
	double phi_min;
	double phi_max;

	// phi_farey is the farey number for this tree node; it is our
	// best guess for the angle at which the zero will occur.
	mpf_t  phi_farey;

	// The actual location of the zero.
	cpx_t z_zero;

	// The actual angle at which the zero is found.
	mpf_t  phi_zero;
	mpf_t  r_zero;

	// The scale difference between phi_farey and phi_zero
	mpf_t phi_offset;

} ZeroNode;

// Global pi;
static mpf_t pi;

void make_zero(int idx, ZeroNode* zn, int ndigits, int nprec)
{
	printf("\n");
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

	// Iniitialize other values;
	cpx_init(zn->z_zero);
	mpf_init(zn->phi_zero);
	mpf_init(zn->r_zero);
	mpf_init(zn->phi_offset);

	// set up guesses for the zero-finder
	double phig = mpf_get_d(zn->phi_farey);
	double phid = 0.5 * (zn->phi_max - zn->phi_min);

	// try to find the zero.
	double rdelta = 0.15;
	int lvl = zn->tree_node.level;
	double r = 1<<(lvl-1);
	if (5 == idx) r = 8.55;
	if (9 == idx) r = 17.4;
	if (10 == idx) r = 22.4;
	if (16 == idx) r = 10.36;
	if (17 == idx) r = 29.19;
	if (23 == idx) r = 24.65;
	if (32 == idx) r = 13.81;
	if (64 == idx) r = 17.759;
	if (128 == idx) r = 22.207;
	for ( ; r< 30; r+= rdelta)
	{
		bool fnd = survey_cell(parti, zn->z_zero, r, phig, rdelta, phid, ndigits, nprec);
		if (fnd)
		{
			bool brk = false;
			mpf_t phi;
			mpf_init(phi);
			fp_arctan2(phi, zn->z_zero->im, zn->z_zero->re, nprec);
			mpf_div(phi, phi, pi);
			double zphi = mpf_get_d(phi);
			if (zn->phi_min < zphi && zphi < zn->phi_max) brk = true;
break;
			mpf_clear(phi);
			if (brk) break;
		}
	}

	// record the radius and the angle.
	cpx_abs(zn->r_zero, zn->z_zero);
	// fp_prt("r = ", zn->r_zero); printf("\n");

	// Well, the angle divided by pi, so it goes from 0.0 to 1.0
	fp_arctan2(zn->phi_zero, zn->z_zero->im, zn->z_zero->re, nprec);
	mpf_div(zn->phi_zero, zn->phi_zero, pi);

	// Compute the mis-placement.
	// If phi_farey < phi_zero, then divide by the difference between
	// farey and right boundary, else the difference to the left boundary.

	mpf_sub(zn->phi_offset, zn->phi_farey, zn->phi_zero);
	mpf_div_ui(zn->phi_offset, zn->phi_offset, 2); // half-angle
printf("---------------\n");
printf("duude angluler diff/2 =%g\n", mpf_get_d(zn->phi_offset));
	int denden = 1;
	if (0.0 < mpf_get_d(zn->phi_offset))
	{
		denden = zn->tree_node.farey.denom * zn->tree_node.left.denom;
	}
	else
	{
		denden = zn->tree_node.farey.denom * zn->tree_node.right.denom;
	}

	mpf_mul_ui(zn->phi_offset, zn->phi_offset, denden);
printf("denden=%d  ", denden);
	fp_prt("offset = ", zn->phi_offset); printf("\n");

	// fp_prt("theta/pi = ", zn->phi_zero); printf("\n");

// #define CHECK_RESULT 1
#ifdef CHECK_RESULT
	cpx_t check;
	cpx_init(check);
	parti(check, zn->z_zero, nprec);
	double eps_r = cpx_get_re(check);
	double eps_i = cpx_get_im(check);
	double eps = sqrt(eps_r*eps_r + eps_i*eps_i);
	printf("fun at zero = %g\n", eps);
	cpx_clear(check);
#endif

	fflush(stdout);

}

void clear_zero(ZeroNode* zn)
{
	mpf_clear(zn->phi_farey);
	cpx_clear(zn->z_zero);
	mpf_clear(zn->phi_zero);
	mpf_clear(zn->r_zero);
	mpf_clear(zn->phi_offset);
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
	int ndigits = 40;
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;
	mpf_set_default_prec(bits);

	// Global initialization
	mpf_init(pi);
	fp_pi(pi, nprec);

	int lvl = atoi(argv[1]);
	int midx = 1<<lvl;

	for (int idx=0; idx < midx; idx++)
	{
		// Skip the right side of the tree.
		int tb = idx;
		while (3 <= tb) { tb >>= 1; }
		if (2 > tb) continue;

		ZeroNode zero;
		make_zero(idx, &zero, ndigits, nprec);

		TreeNode node;
		// make_node(idx, &node);
		node = zero.tree_node;

		printf("node ");
		printf("idx=%d ", idx);
		// printf("lv=(%d, %d) ", node.level, node.vindex);
		printf("dy=%ld/%ld ", node.dyadic.numer, node.dyadic.denom);
		// printf("min=%ld/%ld ", node.left.numer, node.left.denom);
		// printf("max=%ld/%ld ", node.right.numer, node.right.denom);
		printf("far=%ld/%ld ", node.farey.numer, node.farey.denom);
		// printf("phi = (%g %g) ", zero.phi_min, zero.phi_max);
		printf("phx=%g ", mpf_get_d(zero.phi_farey));
		printf("phi=%g ", mpf_get_d(zero.phi_zero));
		printf("r=%g ", mpf_get_d(zero.r_zero));
		printf("off=%g ", mpf_get_d(zero.phi_offset));
		printf("\n");

		clear_zero(&zero);
	}

	printf("Done!\n");
}
