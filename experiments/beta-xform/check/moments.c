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

double coeff[20][20];

double orthonormo(int n, double x)
{
	double sum = 0.0;
	for (int j=0; j<n; j++)
		sum += coeff[n][j] * orthonormo (j, x);
	sum += coeff[n][n] * pow (x, n);
	return sum;
}

double prod(double beta, int n, int m)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double pa = orthonormo(n, tk);
		double pb = orthonormo(m, tk);
		sum += pa * pb / bk;
		bk *= beta;
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
		if (1.0e17 < bk) break;
	}
	return sum;
}

double ortho(double beta, int n, int j)
{
	double sum = 0.0;
	double bk = 1.0;
	double tk = 1.0;
	for (int k=0; k< 453000; k++)
	{
		double tkn = pow(tk, n);
		double ort = orthonormo(j, tk);
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
	double msq = fmoment(beta, 0);
	coeff[0][0] = 1.0 / sqrt(msq);

	for (int n=1; n< max; n++)
	{
		for (int j=0; j<n; j++)
			coeff[n][j] = -ortho(beta, n, j);
		coeff[n][n] = 1.0;
		double msq = prod(beta, n, n);
		double rms = 1.0 / sqrt(msq);
		for (int j=0; j<=n; j++)
			coeff[n][j] *= rms;
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
		int nset = 4;
		setup(beta, nset);
for (int n=0; n<nset; n++) {
for (int m=0; m<=n; m++) {
printf("%d %d co=%f  rat=%f\n", n, m, coeff[n][m], coeff[n][m]/coeff[n][n]);
}}

exit(1);
		printf("%d	%f", i, beta);
		for (int n=0; n<8; n++)
		{
			// double fmom = fmoment(beta, n);
			double fmom = orthonormo(beta, n);
			printf("	%.10f", fmom);
		}
		printf("\n");
		fflush(stdout);
	}
}
