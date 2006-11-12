/* 
 * mp-complex.h
 *
 * Utility routines for handling complex values
 *
 * Linas Vepstas October 2006
 */

#include <gmp.h>

#ifndef __MP_COMPLEX_H__
#define __MP_COMPLEX_H__

typedef struct {
	mpf_t re;
	mpf_t im;
} cpx_t;

static inline void cpx_init (cpx_t *z)
{
	mpf_init (z->re);
	mpf_init (z->im);
}

static inline void cpx_clear (cpx_t *z)
{
	mpf_clear (z->re);
	mpf_clear (z->im);
}

static inline void cpx_add (cpx_t *sum, cpx_t *a, cpx_t *b)
{
	mpf_add (sum->re, a->re, b->re);
	mpf_add (sum->im, a->im, b->im);
}

static inline void cpx_sub (cpx_t *dif, cpx_t *a, cpx_t *b)
{
	mpf_sub (dif->re, a->re, b->re);
	mpf_sub (dif->im, a->im, b->im);
}

static inline void cpx_mul (cpx_t *prod, cpx_t *a, cpx_t *b)
{
	mpf_t pre, pim, tmp;
	mpf_init (pre);
	mpf_init (pim);
	mpf_init (tmp);
	
	mpf_mul (tmp, a->im, b->im);
	mpf_mul (pre, a->re, b->re);
	mpf_sub (pre, pre, tmp);
	
	mpf_mul (tmp, a->im, b->re);
	mpf_mul (pim, a->re, b->im);
	mpf_add (pim, pim, tmp);
	
	mpf_set (prod->re, pre);
	mpf_set (prod->im, pim);

	mpf_clear (pre);
	mpf_clear (pim);
	mpf_clear (tmp);
}

static inline void cpx_recip (cpx_t *recip, cpx_t *z)
{
	mpf_t mag,tmp;
	mpf_init (mag);
	mpf_init (tmp);
	mpf_mul (mag, z->re, z->re);
	mpf_mul (tmp, z->im, z->im);
	mpf_add (mag, mag, tmp);
	mpf_ui_div (mag, 1, mag);
	mpf_mul (recip->re, z->re, mag);
	mpf_mul (recip->im, z->im, mag);
	mpf_neg (recip->im, recip->im);
	
	mpf_clear (mag);
	mpf_clear (tmp);
}

#endif /* __MP_COMPLEX_H__ */
