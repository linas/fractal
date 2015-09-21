
/*
 * Stern-Brocot tree, new algo
 * This is the same old question mark, just a different API.
 *
 * XXX This has been copied to stern-brocot.c  This version
 * is dead.
 *
 * Linas Vepstas September 2015
 */

/**
 * Computes the Stern-Brocot tree.
 * It takes as input a dyadic fraction pos / 2^level
 * assuming 0 <= pos <= 2^level.
 * It returns the corresponding p/q value in p and q.
 *
 * In practical terms, given an input x = pos / 2^level
 * this returns the inverse question-mark, given by
 * (as per usual) ?(p/q) = pos / 2^level.
 */
void stern_brocot_tree(unsigned __int128 pos, int level,
                       unsigned __int128 &p, unsigned __int128 &q);

/** print a 128-bit val; since glibc does not support it. */
int pr128(unsigned __int128);