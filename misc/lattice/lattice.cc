
/*
 * lattice.cc
 *
 * Graph every line of the form (ax+b)/(cx+d) for SL(2,Z)
 * Do the same for the dyadic lines, as well.
 *
 * Linas Sept 2015
 */

struct Uni
{
	// m and k define a dyadic rational m/2^k
	unsigned long m;
	short k;

	// a,b,c,d define a (projective) unimodular matrix.
	long a,b,c,d;

	// multiply into this matrix
	void mult(const Uni&);
	void make_uni(unsigned long, short);
};

void Uni::mult(const Uni& other)
{

}

void Uni::make_uni(unsigned long em, short kay)
{
}

int main(int argc, char* argv[])
{
}
