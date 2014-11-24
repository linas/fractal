//
// schroed.h
//
// Solve schroedinger eqn in a 1D lattice.
//
// Linas Vepstas 23 Nov 2014
//

extern long double *wavefn;

void init(size_t len, long double omega);
long double solve(size_t len, long double omega, long double energy);

void set_pot(size_t len, long double *p);

