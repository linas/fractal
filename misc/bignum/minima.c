
/*
 * minima.c
 *
 * Find the fit for the amplitude of a(n)
 *
 * Linas August 2005
 */


#include <gsl/gsl_multimin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double data [4000];
int num_data_pts = 0;

void
read_data (void)
{
	int c;
	char buf[4000];

	/* read in floating point values in the first column */
	int disc = 1;	
	int i =0;
	int n=0;
	while( (c=getchar()) != EOF)
	{
		if (c == '#')
		{
			printf ("%c", c);
			while( (c=getchar()) != '\n')
			{
				printf ("%c", c);
			}
			printf("\n");
			continue;
		}
		if (disc && c != '\t') continue;
		disc = 0;
		if (c == '\t') continue;

		buf[i] = c;

		if ( c == '\n')
		{
			buf[i] = 0;
			disc = 1;
		   i=-1;	

			double var = atof(buf);

			double corr = exp (-4.0* (sqrt(n+1)));
			var *= corr;
			var /= sin(M_PI*(-2.125+sqrt(2.125*2.125+4.0*(n-1.97)/M_PI)));

// printf ("duude %d %g\n", n, var);
			var = -log (var);
			data[n] = var;
			n++;
		}
		i++;
	}
	num_data_pts = n;
	printf ("read %d datapts\n", num_data_pts);
}

double my_func (double a, double b, double c, double d)
{
	int i;
	double ms = 0.0;
	
	// num_data_pts = 200;
	
	int n=0;
	for (i=7; i<num_data_pts; i++)
	{
		double fitter = a + c*sqrt (1.0+i);
		// double fitter = a + sqrt (b+c*i);
		// double fitter = a + pow (b+c*i, d);
		
		double term = fitter - data[i];
		ms += term*term;
		// printf ("n=%d  %g %g %g %g \n", i, data[i], fitter, term, ms); 
		n++;
	}
	ms /= n;
	// printf ("fit %g %g %g %g gets %g\n", a,b,c, d, ms);
	return ms;
}

double fitter (const gsl_vector * x, void * params)
{
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	double c = gsl_vector_get(x, 2);
	double d = gsl_vector_get(x, 3);

	return my_func (a,b,c, d);
}

void fit (void)
{
	gsl_multimin_fminimizer *fm;

#define FITS 4
	fm = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex, FITS);
	
	gsl_multimin_function f;
	f.f = fitter;
	f.n = FITS;
	f.params = 0;
	
	gsl_vector *start_pt = gsl_vector_alloc (FITS);
	gsl_vector_set (start_pt, 0, 3.0);
	gsl_vector_set (start_pt, 1, 21.0);
	// gsl_vector_set (start_pt, 2, 13.3);
	gsl_vector_set (start_pt, 2, 3.6);
	gsl_vector_set (start_pt, 3, 0.5);
	
	gsl_vector *stepsize = gsl_vector_alloc (FITS);
	gsl_vector_set (stepsize, 0, 0.01);
	gsl_vector_set (stepsize, 1, 0.01);
	gsl_vector_set (stepsize, 2, 0.01);
	gsl_vector_set (stepsize, 3, 0.01);
	
	gsl_multimin_fminimizer_set (fm, &f, start_pt, stepsize); 

	int iter = 0;
	int status;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(fm);
					  
		if (status) 
			break;

		double size = gsl_multimin_fminimizer_size (fm);
		status = gsl_multimin_test_size (size, 1e-5);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}
		printf ("%5d ", iter);
		int i;
		for (i = 0; i <FITS; i++)
		{
			printf ("%11.8g ", gsl_vector_get (fm->x, i));
		}
		printf ("f() = %8.8g size = %.3f\n", fm->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100000);
}

int main (int argc,  char *argv[]) 
{
	read_data ();
	fit();
}

