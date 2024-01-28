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
double proj[20];

double orthonormo(double beta, int n, double x)
{
	if (0 == n) return norm[n];

	double sum = norm[n] * pow(x, n);
	sum += proj[n] * orthonormo(beta, n-1, x);
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

void setup(double beta, int max)
{
	norm[0] = 1.0;
	double msq = prod(beta, 0, 0);
	norm[0] = 1.0 / sqrt(msq);

	for (int n=1; n< max; n++)
	{
		proj[n] = -ortho(beta, n);
		norm[n] = 1.0;
		double msq = prod(beta, n, n);
		double rms = sqrt(msq);
		norm[n] /= rms;
		proj[n] /= rms;
	}
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
		setup(beta, 10);
for (int n=0; n<8; n++) {
for (int m=0; m<=n; m++) {
double n1 = prod(beta, m, n);
printf("%d %d pr=%f\n", m,n,n1);
}}

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
