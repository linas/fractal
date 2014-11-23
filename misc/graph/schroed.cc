//
// schroed.cc
//
// Solve schroedinger eqn ion a 1D lattice.
//
// Linas Vepstas 23 Nov 2014
//

#define LEN 1000
double wavefn[LEN];

double pot[LEN];

void init()
{
	int i;
	for (i=0; i<LEN; i++) wavefn[i] = 0.0;
}


main()
{

	init();
}
