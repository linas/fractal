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
 * Linas Vepstas April 2016
 */

