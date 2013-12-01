
#include <complex.h>
#include <math.h>
#include <stdio.h>

int
main(int argc, char*argv[])
{
	complex lambda=0.2+I*0.3;

	double re = atof(argv[1]);
	double im = atof(argv[2]);

	lambda = re + I* im;

	complex a0 = 1.0;
	complex a1 = -lambda;

	complex anm1 = a1;
	complex anm2 = a0;

	for (int n=2; n<28; n++)
	{
		complex an = (1.0 / ((complex) n)) * (anm2 - lambda*anm1);

		printf("its %d  %g + I %g\n", n, creal(an), cimag(an));
		anm2 = anm1;
		anm1 = an;
	}
}
