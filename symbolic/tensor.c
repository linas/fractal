
/* tensor.c 
 *
 * bogus tensor library because I'm too lazy to learn other systems
 *
 */

#include <glib.h>

/* partial eff /partial x or y */
struct eff_s
{
	int dx;
	int dy;
};

typedef struct eff_s Part;

/* one term */
struct term_s
{
	int coeff;
	GList * parts;   /* list of effs */
};

typedef struct term_s Term;

/* list of terms */
struct expr_s
{
	GList * terms;
};

typedef struct expr_s Expr;

/* =============================================== */

void 
PrintPart (Part *p)
{
	int i;
	if (p->dx != 0 || p->dy != 0) 
	{
		printf ("d");
	}
	for (i=0; i<p->dx; i++)
	{
		printf ("x");
	}
	for (i=0; i<p->dy; i++)
	{
		printf ("y");
	}
	printf ("f ");
}

void
PrintTerm (Term *t)
{
	printf ("+ %d ", t->coeff);
	
	GList *pn;
	for (pn = t->parts; pn; pn=pn->next)
	{
		Part *p = pn->data;
		PrintPart(p);
	}
	printf ("\n");
}

void 
PrintExpr (Expr *e)
{
	printf ("Print expr: \n");

	GList *tn;
	for (tn = e->terms; tn; tn=tn->next)
	{
		Term *t = tn->data;
		PrintTerm (t);
	}

	printf ("The End\n");
}

/* ======================================================== */

Term *Mult_T_T (Term *t1, Term *t2)
{
	Term * t3 = g_new0(Term,1);
	
	t3->coeff = t1->coeff * t2->coeff;
	t3->parts = g_list_copy (t1->parts);
	GList *p = g_list_copy (t2->parts);
	t3->parts = g_list_append (t3->parts, p);

	return t3;
}

/* ======================================================== */

Term *dxf (void)
{
	Term *t = g_new (Term,1);

	Part *p = g_new0 (Part,1);
	p->dx = 1;
	p->dy = 0;

	t->parts = g_list_append (t->parts, p);
	t->coeff = 1;

	return t;
}

Term *dyf (void)
{
	Term *t = g_new (Term,1);

	Part *p = g_new0 (Part,1);
	p->dx = 0;
	p->dy = 1;

	t->parts = g_list_append (t->parts, p);
	t->coeff = 1;

	return t;
}

Expr * setup (void)
{
	Expr *e = g_new0 (Expr,1);

	Term * t = dxf();
	e->terms = g_list_append (e->terms, t);

	return e;
}

main ()
{

	Expr *e = setup();
	PrintExpr (e);
}
