
/*
 * Build a takaagi-style distribution
 *
 * Linas Vepstas October 2008
 */

double triangle (double x)
{
	x -= (int) x;
	if (x < 0.5) return x;
	return 1.0-x;
}

double prod (double x)
{
}

