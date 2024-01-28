/*
 * moments.c
 *
 * Moments.
 * January 2024
 */

#include <math.h> 
#include <stdio.h>
#include <stdlib.h> 

// Hausdorff moment (Hamburger moment on unit interval)
double fmoment(double beta, int n)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow(tk, n);
		sum += tkn / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

double norm[20];
double orth[20];

double orthonormo(double beta, int n, double x)
{
	if (0 == n) return norm[n];

	double sum = norm[n] * pow(x, n);
	sum += orth[n] * orthonormo(beta, n-1, x);
	return sum;
}

double prod(double beta, int n, int m)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double pa = orthonormo(beta, n, tk);
		double pb = orthonormo(beta, m, tk);
		sum += pa * pb / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

double ortho(double beta, int n)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow(tk, n);
		double ort = orthonormo(beta, n-1, tk);
		sum += ort * tkn / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}


int main(int argc, char* argv[])
{
	int imax = 1200;
imax=1;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double beta = x + 1.0;
beta=1.6;

double f0 = fmoment(beta, 0);
printf("f0=%f\n", f0);
norm[0] = 1.0 / sqrt(f0);
double n0 = prod(beta, 0, 0);
printf("n0=%f\n", n0);

double f1 = ortho(beta, 1);
printf("f1=%f\n", f1);

orth[1] = -f1;
norm[1] = 1.0;
double p01 = prod(beta, 0, 1);
printf("p0 x p1=%f\n", p01);

double n1 = prod(beta, 1, 1);
printf("n1=%f\n", n1);
orth[1] /= sqrt(n1);
norm[1] /= sqrt(n1);

double nn1 = prod(beta, 1, 1);
printf("n1=%f\n", nn1);
exit(1);
		printf("%d	%f", i, beta);
		for (int n=0; n<8; n++)
		{
			// double fmom = fmoment(beta, n);
			double fmom = ortho(beta, n);
			printf("	%.10f", fmom);
		}
		printf("\n");
		fflush(stdout);
	}
}
