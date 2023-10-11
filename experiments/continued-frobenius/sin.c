
#include <math.h>
main()
{
	double x;

	for (x=3.141; x<3.142; x+=0.00001)
	{
		double c = 1.5+cos(3.0*x)+0.5*cos(5.0*x);
		double s = sin(3.0*x)+0.5*sin(5.0*x);
		printf ("%f	%g	%g\n", x,c,s);
	}

}
