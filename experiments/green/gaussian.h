/*
 * Semi-Symbolic polynomial handling in c language
 * Very simple, only handles polynomials in one variable
 *
 * Linas Vepstas 16 December 2005
 */

#ifndef GAUSSIAN_H__
#define GAUSSIAN_H__

#define MAXORDER 100

typedef double Poly[MAXORDER];

#define MAXTERMS 100

/* derivatives of the standard normal distribution 
 * with mean mu and standard deviation sigma 
 */
typedef struct 
{
	double sigma;
	double mean;
	double norm;
	Poly derivs[MAXTERMS];
} Gaussian;

Gaussian * gaussian_new (double mu, double sigma);

/* return the 'order'th derivative of a gaussian, at position val */
double gaussian_eval (Gaussian *g, int order, double val);

#endif
