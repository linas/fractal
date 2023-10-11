
double
cube (double u) 
{
	// return u*u*u-2.0*u*u-2.0*u +1.0;
	return 1.0+u*(-2.0+u*(u-2.0));
}

main() {

	int i;
	double ua = 2.0;   // is neg
	double ub = 3.0;   // is pos
	for (i=0; i<40; i++)
	{
		double uc = 0.5 * (ua+ub);
		double fuc = cube (uc);
		if (0.0 > fuc) 
		{
			ua = uc;
		}
		else
		{
			ub = uc;
		}

		printf ("%d: %18.16g %18.16g\n", i,uc,fuc);
	}
}
