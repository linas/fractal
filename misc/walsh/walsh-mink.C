/*
 * walsh-mink.C
 *
 * Dereivative of minkowski-question-mark, integrated with 
 * the Walsh functions, mapped to unit interval.
 *
 * Linas Vepstas August 2009
 */

/*
 * argument is interpreted as p/2^n
 */
double eval (int p, int n)
{
	return 0.0;
}

int main(int argc, char * argv[])
{

	int n = 4;
	int m = 1<<n:
	for (int p=0; p<m; p++)
	{
		double x = ((double) p) / ((double) m);
		double y = eval(q, n);
		printf("%d	%f	%f\n", p, x, y);
	} 
}
