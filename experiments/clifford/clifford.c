/**
 * clifford.c
 *
 * Finally getting around to exploring a very old idea:
 * The basis of a clifford algebra is in teh form of a binary tree.
 * Given that this project is an exploration of binary trees, we
 * may as well take a look at what happens here. We take the vector
 * space V as being infinite-dimansional, so that the tree is
 * infinite-dimensional, too.
 *
 * We have a (non-canonical) mapping of the clifford algebra into a
 * tree.  This mapping creates a bijection between the clifford algebra,
 * and functions on the dyadic rationals.  Multiplication in the
 * clifford algebra thereby induces multiplication of such trees.
 * The goal here is to explore "what this looks like".
 *
 * The mapping is the "bitstring mapping".  If b_0 b_1 b_2... is a
 * bitstring, then the correpsonding clifford-algebra basis function
 * is given by the indexes of those bits that are not zero.  That is,
 * if the non-zero bits are b_i b_j b_k ... then the corresponding
 * clifford alg element is e_i e_j e_k ...
 *
 * Unfinished....
 * Clifford algebras are quantizations of exterior algebras and are
 * therefore kind-of-ish homological in thier behavior.  The homological
 * variant is explored in the directory `../alternate`.
 *
 * Unfinished because it requires a tedious implementation of a product,
 * which is hard to do in C... because of the lack of lambdas in C.
 *
 * Linas Vepstas April 2016
 */


#include <math.h>
#include <stdio.h>

#define BITLEN 56
typedef char bitstring[BITLEN];

void to_bits(bitstring bits, double x)
{
	x -= floor(x);

	for (int i=0; i<BITLEN; i++)
	{
		if (0.5 < x) { bits[i] = 1; x -= 0.5; }
		else bits[i] = 0;
		x *= 2.0;
	}
}

double from_bits(bitstring bits)
{
	double val = 0.0;
	double bv = 0.5;
	for (int i=0; i<BITLEN; i++)
	{
		if (bits[i]) val += bv;
		bv *= 0.5;
	}
	return val;
}

/* An element of the clifford algebra is a vector i.e. a function
 * f that takes on a real value for any given basis vector in the
 * cliff alg, and since any given basis vector is represented by
 * a bit-string, that's what we do here.
 */
typedef double (*cliff)(bitstring);

/* product of two vectors */
double product(cliff f, cliff g, bitstring b)
{


}


int main(int argc, char* argv[])
{

}
