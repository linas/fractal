
/*
 * cont_frac.c
 *
 * Convert real numbers to continued fractions, and print them.
 * Uses arbitrary-precision math to do so.
 *
 * Linas Vepstas April 2006
 *
 */

#include <gmp.h>
#include <malloc.h>
#include <stdio.h>

typedef struct 
{
	int nterms;
	int len;
	unsigned long *tinued_frac;
} continued_frac_t;


continued_frac_t *continued_frac_new (int len)
{
	continued_frac_t *cf;
	cf = (continued_frac_t *) malloc (sizeof (continued_frac_t));
	cf->len = len;
	cf->nterms = 0;
	cf->tinued_frac = (unsigned long *) malloc (len*sizeof(unsigned long));
	return cf;
}

// convert mp float to a continued fraction
void continued_frac_from_mpf (continued_frac_t *cf, mpf_t x)
{
	int i;
	mpf_t val, inv;
	mpf_init (val);
	mpf_init (inv);

	mpf_set (val, x);
	for (i=0; i<cf->len; i++)
	{
		unsigned long int_part = mpf_get_ui (val);
		cf->tinued_frac[i] = int_part;

		if (0 == int_part && i) break;
		mpf_sub_ui (val, val, int_part);
		mpf_ui_div (inv, 1, val);
		mpf_set (val, inv);
	}
	cf->nterms = i;

	mpf_clear (val);
	mpf_clear (inv);
}

// convert double to continued fraction
void continued_frac_from_double (continued_frac_t *cf, double x)
{
	mpf_t val;
	mpf_init (val);
	mpf_set_d (val, x);
	continued_frac_from_mpf (cf, val);
	mpf_clear (val);
}

// convert string to continued fraction
void continued_frac_from_str (continued_frac_t *cf, char *str)
{
	mpf_t val;
	mpf_init (val);
	mpf_set_str (val, str, 10);
	continued_frac_from_mpf (cf, val);
	mpf_clear (val);
}

// convert continued fraction to mpf
void continued_frac_to_mpf (mpf_t result, continued_frac_t *cf)
{
	mpf_t one;

	mpf_init (one);
	mpf_set_ui (one, 1);
	int i;
	mpf_set_ui (result, 0);
	for (i=cf->nterms-1; i>0; i--)
	{
		mpf_add_ui (result, result, cf->tinued_frac[i]);
		mpf_div (result, one, result);
	}
	mpf_add_ui (result, result, cf->tinued_frac[0]);

	mpf_clear (one);
}

// print terms of continued fraction
void continued_frac_print (continued_frac_t *cf)
{
	int i;
	printf ("nterms = %d\n", cf->nterms);
	for (i=0; i<cf->nterms; i++)
	{
		printf ("its %d %lu \n", i, cf->tinued_frac[i]);
	}
}

void fp_prt (char * str, mpf_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 60, val);
	printf ("\n");
}


int main ()
{
	continued_frac_t *cf = continued_frac_new(100);

#if 0
	continued_frac_from_double (cf, 3.14159265358979);
	continued_frac_from_str (cf, "3.14159265358979");
	continued_frac_print (cf);

	mpf_t pi;
	mpf_init (pi);
	
	int k;
	for (k=1; k<30; k++)
	{
		cf->nterms = k;
		continued_frac_to_mpf (pi, cf);
		fp_prt (" yo " , pi);
	}
#endif 

#if 0
	continued_frac_from_str (cf, "0.772156649015328606065120900824024310421593359399235988057672e-1");
	continued_frac_print (cf);

	continued_frac_from_str (cf, "0.474863147741964237030450675938983643268438533652437993798812e-2");
	continued_frac_print (cf);

	continued_frac_from_str (cf, "0.366100893495480732657022726472761255128467820041878332367143e-3");
	continued_frac_print (cf);
	
	continued_frac_from_str (cf, "0.376007305663723248040089177144004381633628176123569592709274e-3");
	continued_frac_print (cf);
#endif
	
	continued_frac_from_str (cf, "0.14301182486231445523918610317570169576931647768652503561046e-3");
	continued_frac_print (cf);
	
#if 0
	continued_frac_from_str (cf, "
	continued_frac_print (cf);
#endif
	
	return 0;
}
