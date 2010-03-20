
/*
 * Attempt a linear decomposition of continuous/differentiable
 * eigenvalues into fractals.
 *
 * March 2010
 */

double triangle(double x)
{
	x -= floor(x);
	x *= 2.0;
	if (x < 1.0) return x;
	return 2.0-x;
}

double tagaki(double w, double x)
{
	int i;
	double tn = 1.0
	double wn = 1.0;
	double acc = 0.0;

	for (i=0; i< 50; i++)
	{
		acc += wn * triangle(tn*x);
		tn *= 2.0;
		wn *= w;
	}

	return acc;
}
