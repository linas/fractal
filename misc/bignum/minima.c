
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

			var = -log (var);
			data[n] = var;
			n++;
		}
		i++;
	}
	num_data_pts = n;
}

double my_func (double a, double b, double c)
{
	int i;
	double ms = 0.0;
	
	for (i=0; i<num_data_pts; i++)
	{
		fitter = a + sqrt (b+c*i);
		
		double term = fitter - data[i];
		ms += term*term;
	}
	return ms;
}

double fitter (const gsl_vector * x, void * params)
{
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	double c = gsl_vector_get(x, 2);

	return my_func (a,b,c);
}

void fit (void)
{
	gsl_multimin_fminimizer *fm;

	fm = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex, 3);
	
	gsl_multimin_function f;
	f.f = fitter;
	f.n = 3;
	f.params = 0;
	
	gsl_vector *start_pt = gsl_vector_alloc (3);
	gsl_vector_set (start_pt, 0, 5.0);
	gsl_vector_set (start_pt, 1, 5.0);
	gsl_vector_set (start_pt, 2, 5.0);
	
	gsl_vector *stepsize = gsl_vector_alloc (3);
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

		double size = gsl_multimin_fminimizer_size (fm);
		status = gsl_multimin_test_size (size, 1e-4);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}
		printf ("%5d ", iter);
		int i;
		for (i = 0; i <3; i++)
		{
			printf ("%10.3e ", gsl_vector_get (fm->x, i));
		}
		printf ("f() = %7.3f size = %.3f\n", fm->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);
}

int main (int argc,  char *argv[]) 
{
	fit();
}

