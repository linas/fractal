
/*
 * minima.c
 *
 * Find the fit for the amplitude of a(n)
 *
 * Linas August 2005
 */


#include <gsl/gsl_multimin.h>
#include <math.h>

double my_func (double a, double b, double c)
{
	return (a-1)*(a-1) + (b-2)(b-2) + (c-3)*(c-3);
}

double fitter (const gsl_vector * x, void * params)
{
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	double c = gsl_vector_get(x, 2);

	my_func (a,b,c);
}

void fit (void)
{
	gsl_multimin_fminimizer *fm;

	fm = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex, 3);
	
	gsl_multimin_function f;
	f.f = fitter;
	f.n = 3;
	f.params = 0;
	
	gsl_vector start_pt = gsl_vector_alloc (3);
	gsl_vector_set (start_pt, 0, 5.0);
	gsl_vector_set (start_pt, 1, 5.0);
	gsl_vector_set (start_pt, 2, 5.0);
	
	gsl_vector stepsize = gsl_vector_alloc (3);
	gsl_vector_set (stepsize, 0, 0.01);
	gsl_vector_set (stepsize, 1, 0.01);
	gsl_vector_set (stepsize, 2, 0.01);
	
	gsl_multimin_fminimizer_set (fm, &f, start_pt, stepsize); 

	int iter = 0;
	int status;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(fm);
					  
		if (status) 
			break;

		int size = gsl_multimin_fminimizer_size (fm);
		status = gsl_multimin_test_size (size, 1e-2);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}
		printf ("%5d ", iter);
		for (i = 0; i <3; i++)
		{
			printf ("%10.3e ", gsl_vector_get (fm->x, i));
		}
		printf ("f() = %7.3f size = %.3f\n", fm->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);
}

main () 
{
	fit();
}
