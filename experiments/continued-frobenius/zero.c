

/*
 * zero.c
 *
 * zero games.
 */


#include <complex.h>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// return binomial coefficient s! / (n! (s-n)!)
double complex cbin (double complex s, int n)
{
	int k;

	long double complex bin = 1.0L;
	for (k=1; k<=n; k++)
	{
		bin *= s / ((long double) k);
		s -= 1.0L;
	}
	return bin;
}


#define NELT 8

void init_arr (double *arr)
{
	arr[0] = 1.0;
	arr[1] = 0.5;
	arr[2] = 0.4;
	arr[3] = 0.3;
	arr[4] = 0.2345678;
	arr[5] = 0.1987654;
	arr[6] = 0.13;
	arr[7] = 0.114;
}

double complex 
fun (double complex ess, double *arr)
{
	double complex val = 0.0;

	int i;
	for (i=0; i<NELT; i++)
	{
		val += cbin (ess, i) * arr[i];
	}
	return val;
}

int
my_what_fun (const gsl_vector * x, void * p, gsl_vector * f) 
{
	double *arr = (double *) p;
	const double re = gsl_vector_get(x,0);
	const double im = gsl_vector_get(x,1);

	const double complex ess = re + I*im;
	double complex val = fun (ess, arr);

	gsl_vector_set (f, 0, creal (val));
	gsl_vector_set (f, 1, cimag (val));
	return GSL_SUCCESS;
}

void doit(double * arr)
{
	const gsl_multiroot_fsolver_type * solver_type 
			      = gsl_multiroot_fsolver_hybrid;
// 			      = gsl_multiroot_fsolver_hybrids;

	/* we want a 2d solver */
	gsl_multiroot_fsolver * solver 
			      = gsl_multiroot_fsolver_alloc (solver_type, 2);

	// function declaration
	gsl_multiroot_function gmf = {my_what_fun, 2, (void *)arr};

	// initial guess 
	gsl_vector *guess = gsl_vector_alloc (2);

#define DOM 30.0
	double re, im;
	for (re = -DOM; re < DOM; re +=1.0)
	{
		for (im = 0.0; im < DOM; im +=1.0)
		{
			gsl_vector_set (guess, 0, re);
			gsl_vector_set (guess, 1, im);

			gsl_multiroot_fsolver_set (solver, &gmf, guess);
	
			int i;
			for (i=0; i<100; i++)
			{
				int status = gsl_multiroot_fsolver_iterate (solver);
		
				if (status) break; // solver is stuck
				
				status = gsl_multiroot_test_residual (solver->f, 1e-10);
				if (GSL_CONTINUE != status) break;
			}
	
			gsl_vector * root = gsl_multiroot_fsolver_root (solver);
	
			double rre = gsl_vector_get(root,0);
			double rim = gsl_vector_get(root,1);
	
			// if (fabs (im) < 1e-8) im = 0.0;
			printf ("duude re=%f im=%f\n", rre, rim);
		}
	}
}



int
main (int argc, char *argv[])
{
	double arr[NELT];

	init_arr(arr);

	arr[0] = atof (argv[1]);
	printf ("trying %f\n", arr[0]);

	doit (arr);
	return 0;
}
