/*
 * present.c
 *
 * hunt out group presentations
 *
 * Linas February 2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _matrix {
	int dim;
	int *elts;
} Matrix;

Matrix * matrix_new (int dim)
{
	Matrix *mat = (Matrix *)malloc (sizeof (Matrix));
	mat->dim = dim;
	mat->elts = (int *)malloc (dim*dim*sizeof(int));
	return mat;
}

/* element acccessor */
#define MELT(mat,m,n) ((mat)->elts[(m)*(mat)->dim + (n)])

void matrix_set (Matrix *mat, int m, int n, int val)
{
	MELT(mat,m,n) = val;
}

void matrix_copy (Matrix *to, Matrix *from)
{
	int m,n;
	for (m=0; m<from->dim; m++)
	{
		for (n=0; n<from->dim; n++)
		{
			MELT(to,m,n) = MELT(from, m,n);
		}
	}
}

Matrix * matrix_dup (Matrix *from)
{
	Matrix *mat = matrix_new (from->dim);
	matrix_copy (mat, from);
	return mat;
}

void matrix_unit (Matrix *mat)
{
	int m,n;
	for (m=0; m<mat->dim; m++)
	{
		for (n=0; n<mat->dim; n++)
		{
			MELT(mat,m,n) = 0;
		}
		MELT(mat,m,m) = 1;
	}
}

void matrix_mult (Matrix *prod, Matrix *a, Matrix *b)
{
	int m,n,p;
	for (m=0; m<a->dim; m++)
	{
		for (n=0; n<a->dim; n++)
		{
			int acc = 0;
			for (p=0; p<a->dim; p++)
			{
				acc += MELT(a,m,p) * MELT(b, p,n);
			}
			MELT(prod,m,n) = acc;
		}
	}
}

static inline int matrix_equal (Matrix *a, Matrix *b)
{
	int m,n;
	for (m=0; m<a->dim; m++)
	{
		for (n=0; n<a->dim; n++)
		{
			if (MELT(a,m,n) != MELT(b, m,n)) return 0;
		}
	}
	return 1;
}

static inline int matrix_projective_equal (Matrix *a, Matrix *b)
{
	if (matrix_equal(a,b)) return 1;

	int m,n;
	for (m=0; m<a->dim; m++)
	{
		for (n=0; n<a->dim; n++)
		{
			if (MELT(a,m,n) != -MELT(b, m,n)) return 0;
		}
	}
	return 1;
}

/* ---------------------------------------------------- */
/* linked list of matrixes */
typedef struct _mlist MatList;

struct _mlist {
	MatList *next;  // next in list
	Matrix *mat;    // this group element
	char *str;      // word describing this group element
	char *last;     // last letter of the word
	MatList *match; // equivalent word
};

MatList *matlist_prepend (MatList *ml, Matrix *mat, char *str, char letter)
{
	MatList *node = (MatList *)malloc(sizeof (MatList));
	node->mat = mat;
	node->next = ml;

	/* append letter to string, store string */
	int len = strlen(str);
	node->str = malloc (len+2);
	strcpy(node->str, str);
	node->str[len] = letter;
	node->str[len+1] = 0x0;
	node->last = &node->str[len];

	node->match = NULL;
	return node;
}

/* Find a matching matrix in the list */
MatList *matlist_find (MatList *ptr, Matrix *mat)
{
	while (ptr)
	{
		// if (matrix_equal (ptr->mat, mat)) return ptr;
		if (matrix_projective_equal (ptr->mat, mat)) return ptr;
		ptr = ptr->next;
	}
	return NULL;
}

int matlist_len (MatList *ptr)
{
	int len = 0;
	while (ptr)
	{
		len++;
		ptr = ptr->next;
	}
	return len;
}

/* ---------------------------------------------------- */
/* one-level deep hash table of matrixes */
typedef struct _htab HashTab;

struct _htab {
	int size;
	MatList **matlist;  // hash table
};

HashTab * hashtab_new (int size)
{
	int i;
	HashTab *htab;
	htab = (HashTab *) malloc(sizeof (HashTab));
	htab->size = size;
	htab->matlist = (MatList **) malloc(sizeof (MatList *));
	for (i=0; i<size;i++)
	{
		htab->matlist[i] = NULL;
	}
	return htab;
}

MatList * hashtab_add (HashTab *htab, Matrix *mat, char *str, char letter)
{
	int hash = MELT (mat, 0, 0);
	hash %= htab->size;
	if (0 > hash) hash = -hash;

	MatList *ptr = htab->matlist[hash];
	ptr = matlist_prepend (ptr, mat, str, letter);
	htab->matlist[hash] = ptr;
	return ptr;
}

/* Find a matching matrix in the list */
MatList *hashtab_find (HashTab *htab, Matrix *mat)
{
	int hash = MELT (mat, 0, 0);
	hash %= htab->size;
	if (0 > hash) hash = -hash;

	MatList *ptr = htab->matlist[hash];
	ptr = matlist_find (ptr, mat);
	return ptr;
}

void hashtab_report (HashTab *htab)
{
	int i;
	for (i=0; i<htab->size;i++)
	{
		int len = matlist_len (htab->matlist[i]);
		printf ("bucket %d len %d\n",i,len);
	}
}

/* ---------------------------------------------------- */
/* list of strings */

typedef struct _slist StrList;

struct _slist {
	StrList *next;
	char * str;
	size_t len;
};

StrList * slist_prepend (StrList *sl, char * str)
{
	StrList *head = (StrList *)malloc (sizeof(StrList));
	head->str = str;
	head->len = strlen (str);
	head->next = sl;
	return head;
}

/* Turn a matching pair into a single word string equal to the identity */
char * make_str (MatList *word)
{
	int la = strlen (word->str);
	int lb = strlen (word->match->str);
	lb --; // strip off leading E
	la --;

	char *e = malloc (la+lb+1);
	strcpy (e, &word->str[1]);   // be sure to srip off leading E
	while (lb)
	{
		e[la] = 0x20 ^ word->match->str[lb];
		la++;
		lb--;
	}
	e[la] = 0x0;
	return e;
}

/* Check to see if string a is a subset of string b starting at offset 
 * Cyclic embedding is assumed. Assume that b is equal or longer than a */
static int is_cyclic_off (char *a, char * b, size_t offset)
{
	char * p = b+offset;
	while (*a)
	{
		if (*a != *p) return 0;
		a++;
		p++;
		if (0x0 == *p) p = b;
	}
	return 1;
}

/* Check to see if string a shows up as an inverted substring of string b. 
 * Cyclic embedding is assumed. Assume b is equal or longer than a */
static int is_cyclic_rev (char *a, char * b, size_t offset)
{
	char * p = b+offset;
	while (*a)
	{
		if (*a != (0x20^(*p))) return 0;
		a++;
		if (b == p) p = b+strlen(b);
		p--;
	}
	return 1;
}

/* Check to see if string a or its revers is ebedded in string b.
 * len should the the length of b, which should be equal to or longer than a.
 */
int is_cyclic (char *a, char * b, size_t len)
{
	size_t i;
	for (i=0; i<len; i++)
	{
		if (1 == is_cyclic_off (a,b,i)) return 1;
		if (1 == is_cyclic_rev (a,b,i)) return 1;
	}
	return 0;
}

/* If the word contains a substring that is in the list of known words, 
 * then return true, else false */
int slist_in_list (StrList *sl, char * word)
{
	int len = strlen (word);
	while (sl)
	{
		if (is_cyclic (sl->str, word, len)) return 1;
		sl= sl->next;
	}
	return 0;
}

/* ---------------------------------------------------- */
/* trim -- return new string, equal to a, with letters starting at 
 * off removed, running for length len */
char * trim_off (char * a, size_t off, size_t len)
{
	size_t olen = strlen(a);
	char * cp = malloc (olen-len+1);
	if (off+len <= olen)
	{
		strncpy(cp,a,off);
		cp[off] = 0x0;
		a += off+len;
		strcat (cp, a);
		return cp;
	}
	else
	{
printf ("not implemmented\n");
	}
	return NULL;
}

/* if string a appears inside of string b, then a new string is returned, 
 * where a has been cut out from b. The len argument should be length of b.
 */
char * do_trim (char *a, char * b, size_t len)
{
	size_t i;
	for (i=0; i<len; i++)
	{
		if (1 == is_cyclic_off (a,b,i))
		{
			return trim_off (b, i, strlen (a)); 
		}
	}
	return NULL;
}

/* If the word contains a substring that already appears in the list of known words, 
 * then the word is trimmed. */
char * do_trim_to_slist (StrList *sl, char * word)
{
	int len = strlen (word);
	while (sl)
	{
		if (len > sl->len)
		{
printf ("duude try this %s in %s ??\n", sl->str, word);
			char * shorte = do_trim (sl->str, word, len);
			if (shorte) return shorte;
		}
		sl= sl->next;
	}
	return NULL;
}

char * trim_to_slist (StrList *sl, char * word)
{
	char * shorte = NULL;
	while (1)
	{
		char *even_shorter = do_trim_to_slist(sl,word);
		if (NULL == even_shorter) return shorte;

		if (shorte) free(shorte);
		shorte = even_shorter;
		word = shorte;
	}
	return shorte;
}

/* ---------------------------------------------------- */
/* master struct */

typedef struct _present {
	int dim;
	MatList *generators;    // list of group generators (and thier inverses)
	// MatList *words;         // list of words tested so far 
	HashTab *words;         // hash table of words tested so far 
	MatList *presentation;  // presentation of group so far.
	int found;
	int cnt;
	StrList *unique;
	MatList *ident;         // bogus -- pointer to identity matrix
} Present;

Present * present_new (int dim)
{
	Present *pr = (Present *) malloc (sizeof (Present));
	pr->generators = NULL;
	pr->words = hashtab_new (13);
	pr->presentation = NULL;
	pr->cnt =0;
	pr->found = 0;
	pr->unique = NULL;

	/* identity matrix */
	Matrix *e = matrix_new (dim);
	matrix_unit (e);
	// pr->words = matlist_prepend (NULL, e, "", 'E');
	pr->ident = hashtab_add (pr->words, e, "", 'E');
	pr->dim = dim;

	return pr;
}
	
static inline isinv (char a, char b)
{
	if ((a^0x20) == b) return 1;
	return 0;
}

void present_walk_tree (Present *pr, MatList *node, int depth)
{
	if (0 == depth) return;
	depth --;
	pr->cnt ++;

	MatList *gen = pr->generators;
	while (gen)
	{
		/* Do not concatent matrix to its inverse -- skip these words */
		if (isinv (node->last[0], gen->str[0]))
		{
			gen = gen->next;
			continue;
		}

		/* concatenate matrixes, make a new word */
		Matrix *mat = matrix_new (node->mat->dim);
		matrix_mult (mat, node->mat, gen->mat);

		/* Have we seen this word before ? */
		// MatList *match = matlist_find (pr->words, mat);
		MatList *match = hashtab_find (pr->words, mat);
		if (match)
		{
			pr->presentation = matlist_prepend (pr->presentation, mat, node->str, gen->str[0]);

			/* add forward and backward links */
			match->match = pr->presentation;
			pr->presentation->match = match;
			pr->found ++;

			printf ("got one! %d %s == (%s %c)\n", pr->found, match->str, node->str, gen->str[0]);
		}
		else
		{
			/* unknown. keep going */
			// pr->words = matlist_prepend (pr->words, mat, node->str, gen->str[0]);
			MatList * unk = hashtab_add (pr->words, mat, node->str, gen->str[0]);
			present_walk_tree (pr, unk, depth);
		}

		gen = gen->next;
	}
}

void present_cleanup (Present *pr)
{
	MatList *word = pr->presentation;
	while (word)
	{
		char * sword = make_str (word);
		if (0 == slist_in_list(pr->unique, sword))
		{
printf ("found unique %s\n", sword);
			pr->unique = slist_prepend (pr->unique, sword);
		}
		else
		{
			char * shorte = trim_to_slist (pr->unique, sword);
			if (shorte)
			{
printf ("found short=%s\n", shorte);
			}
			free(sword);
		}

		word = word->next;
	}
}

/* ---------------------------------------------------- */
/* set up braid group 3 */

Present * setup_braid3 (void)
{
	Present *pr = present_new(2);
	
	MatList *ml = NULL;

	/* Set up sigma 1 and sigma 2 of braid group B_3 */
	Matrix *s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = 1;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'A');
	
	Matrix *s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = -1;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'B');
	
	/* and now thier inverses (by hand) */
	s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = -1;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'a');
	
	s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = 1;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'b');

	pr->generators = ml;
	return pr;
}

/* ---------------------------------------------------- */
/* set up free subgroup of braid group 3 */

Present * setup_braid3free (void)
{
	Present *pr = present_new(2);
	
	MatList *ml = NULL;

	/* Set up A = (sigma_1)^2 and B=(sigma_2)^2 of braid group B_3 */
	Matrix *s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = 2;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'A');
	
	Matrix *s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = -2;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'B');
	
	/* and now thier inverses (by hand) */
	s1 = matrix_new (2);
	MELT(s1, 0,0) = 1;
	MELT(s1, 0,1) = -2;
	MELT(s1, 1,0) = 0;
	MELT(s1, 1,1) = 1;
	ml = matlist_prepend (ml, s1, "", 'a');
	
	s2 = matrix_new (2);
	MELT(s2, 0,0) = 1;
	MELT(s2, 0,1) = 0;
	MELT(s2, 1,0) = 2;
	MELT(s2, 1,1) = 1;
	ml = matlist_prepend (ml, s2, "", 'b');

	pr->generators = ml;
	return pr;
}

/* ---------------------------------------------------- */
/* set up 3D Heisenberg group */

Present * setup_heisenberg (void)
{
	Present *pr = present_new(3);
	
	MatList *ml = NULL;

	/* Set up x and y of Heisenberg group */
	Matrix *x = matrix_new (3);
	MELT(x, 0,0) = 1;
	MELT(x, 0,1) = 1;
	MELT(x, 0,2) = 0;
	MELT(x, 1,0) = 0;
	MELT(x, 1,1) = 1;
	MELT(x, 1,2) = 0;
	MELT(x, 2,0) = 0;
	MELT(x, 2,1) = 0;
	MELT(x, 2,2) = 1;
	ml = matlist_prepend (ml, x, "", 'X');

	Matrix *y = matrix_new (3);
	MELT(y, 0,0) = 1;
	MELT(y, 0,1) = 0;
	MELT(y, 0,2) = 0;
	MELT(y, 1,0) = 0;
	MELT(y, 1,1) = 1;
	MELT(y, 1,2) = 1;
	MELT(y, 2,0) = 0;
	MELT(y, 2,1) = 0;
	MELT(y, 2,2) = 1;
	ml = matlist_prepend (ml, y, "", 'Y');

	x = matrix_new (3);
	MELT(x, 0,0) = 1;
	MELT(x, 0,1) = -1;
	MELT(x, 0,2) = 0;
	MELT(x, 1,0) = 0;
	MELT(x, 1,1) = 1;
	MELT(x, 1,2) = 0;
	MELT(x, 2,0) = 0;
	MELT(x, 2,1) = 0;
	MELT(x, 2,2) = 1;
	ml = matlist_prepend (ml, x, "", 'x');

	y = matrix_new (3);
	MELT(y, 0,0) = 1;
	MELT(y, 0,1) = 0;
	MELT(y, 0,2) = 0;
	MELT(y, 1,0) = 0;
	MELT(y, 1,1) = 1;
	MELT(y, 1,2) = -1;
	MELT(y, 2,0) = 0;
	MELT(y, 2,1) = 0;
	MELT(y, 2,2) = 1;
	ml = matlist_prepend (ml, y, "", 'y');

	pr->generators = ml;
	return pr;
}
	
/* ---------------------------------------------------- */
/* set up trilog (polylog n=3) */

Present * setup_trilog (void)
{
	Present *pr = present_new(4);
	
	MatList *ml = NULL;
	Matrix *g0 = matrix_new (4);
	matrix_unit (g0);
	MELT(g0,0,1) = 1;
	MELT(g0,0,2) = 1;
	MELT(g0,1,2) = 2;
	ml = matlist_prepend (ml, g0, "", 'A');
	
	Matrix *g1 = matrix_new (4);
	matrix_unit (g1);
	MELT(g1,2,3) = 1;
	ml = matlist_prepend (ml, g1, "", 'B');
	
	g0 = matrix_new (4);
	matrix_unit (g0);
	MELT(g0,0,1) = -1;
	MELT(g0,0,2) = 1;
	MELT(g0,1,2) = -2;
	ml = matlist_prepend (ml, g0, "", 'a');
	
	g1 = matrix_new (4);
	matrix_unit (g1);
	MELT(g1,2,3) = -1;
	ml = matlist_prepend (ml, g1, "", 'b');
	
	pr->generators = ml;
	return pr;
}

/* ---------------------------------------------------- */
/* set up polylog (polylog n=4) */

Present * setup_quadlog (void)
{
	Present *pr = present_new(5);
	
	MatList *ml = NULL;
	Matrix *g0 = matrix_new (5);
	matrix_unit (g0);
	MELT(g0,0,1) = 1;
	MELT(g0,0,2) = 1;
	MELT(g0,1,2) = 2;
	MELT(g0,0,3) = 1;
	MELT(g0,1,3) = 3;
	MELT(g0,2,3) = 3;
	ml = matlist_prepend (ml, g0, "", 'A');
	
	Matrix *g1 = matrix_new (5);
	matrix_unit (g1);
	MELT(g1,3,4) = 1;
	ml = matlist_prepend (ml, g1, "", 'B');
	
	g0 = matrix_new (5);
	matrix_unit (g0);
	MELT(g0,0,1) = -1;
	MELT(g0,0,2) = 1;
	MELT(g0,1,2) = -2;

	MELT(g0,0,3) = -1;
	MELT(g0,1,3) = 3;
	MELT(g0,2,3) = -3;
	ml = matlist_prepend (ml, g0, "", 'a');
	
	g1 = matrix_new (5);
	matrix_unit (g1);
	MELT(g1,3,4) = -1;
	ml = matlist_prepend (ml, g1, "", 'b');
	
	pr->generators = ml;
	return pr;
}

/* ---------------------------------------------------- */
/* set up polylog (polylog n=5) */

Present * setup_pentalog (void)
{
	Present *pr = present_new(6);

	MatList *ml = NULL;
	Matrix *g0 = matrix_new (6);
	matrix_unit (g0);
	MELT(g0,0,1) = 1;
	MELT(g0,0,2) = 1;
	MELT(g0,0,3) = 1;
	MELT(g0,0,4) = 1;
	MELT(g0,1,2) = 2;
	MELT(g0,1,3) = 3;
	MELT(g0,2,3) = 3;
	MELT(g0,1,4) = 4;
	MELT(g0,2,4) = 6;
	MELT(g0,3,4) = 4;
	ml = matlist_prepend (ml, g0, "", 'A');
	
	Matrix *g1 = matrix_new (6);
	matrix_unit (g1);
	MELT(g1,4,5) = 1;
	ml = matlist_prepend (ml, g1, "", 'B');
	
	g0 = matrix_new (6);
	matrix_unit (g0);
	MELT(g0,0,1) = -1;

	MELT(g0,0,2) = 1;
	MELT(g0,1,2) = -2;

	MELT(g0,0,3) = -1;
	MELT(g0,1,3) = 3;
	MELT(g0,2,3) = -3;

	MELT(g0,0,4) = 1;
	MELT(g0,1,4) = -4;
	MELT(g0,2,4) = 6;
	MELT(g0,3,4) = -4;
	ml = matlist_prepend (ml, g0, "", 'a');
	
	g1 = matrix_new (6);
	matrix_unit (g1);
	MELT(g1,4,5) = -1;
	ml = matlist_prepend (ml, g1, "", 'b');
	
	pr->generators = ml;
	return pr;
}

/* ---------------------------------------------------- */

main ()
{
	int depth;
	Present *pr;

	// pr = setup_braid3();
	pr = setup_braid3free();
	// pr = setup_heisenberg();
	// pr = setup_trilog();
	// pr = setup_quadlog();
	// pr = setup_pentalog();

	for (depth=2; depth <17; depth++)
	{
		present_walk_tree (pr, pr->ident, depth);
		// present_cleanup (pr);

		printf ("at depth %d tested %d words\n", depth, pr->cnt);

		/* now do another go-around */
		hashtab_report (pr->words);
		pr->presentation = NULL;  // XXX should be a free 
		pr->words = NULL; // XXX should be free

		/* identity matrix */
		pr->words = hashtab_new(133);
		Matrix *e = matrix_new (pr->dim);
		matrix_unit (e);
		pr->ident = hashtab_add (pr->words, e, "", 'E');
	}
}
