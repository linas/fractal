
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
	if (0 == t->coeff) return;

	if (1 == t->coeff && NULL != t->parts)
	{
		printf ("+ ");
	}
	else if (0 > t->coeff)
	{
		printf ("- %d ", - t->coeff);
	}
	else 
	{
		printf ("+ %d ", t->coeff);
	}
	
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
	if (NULL == e) 
	{
		printf ("NULL expresion!\n");
		return;
	}

	GList *tn;
	for (tn = e->terms; tn; tn=tn->next)
	{
		Term *t = tn->data;
		PrintTerm (t);
	}

	printf ("The End\n\n");
}

/* ======================================================== */

/* multipliy terms */

Term *Mult_T_T (Term *t1, Term *t2)
{
	Term * t3 = g_new0(Term,1);
	
	t3->coeff = t1->coeff * t2->coeff;
	t3->parts = NULL;

	GList *pn;
	for (pn = t1->parts; pn; pn=pn->next)
	{
		Part *p = pn->data;
		t3->parts = g_list_append (t3->parts, p);
	}
	for (pn = t2->parts; pn; pn=pn->next)
	{
		Part *p = pn->data;
		t3->parts = g_list_append (t3->parts, p);
	}

	return t3;
}

void Mult_T_I (Term *t, int c)
{
	t->coeff *= c;
}

/* ======================================================== */

Part * 
Copy_Part (Part *p)
{
	Part *p3 = g_new0(Part,1);
	p3->dx = p->dx;
	p3->dy = p->dy;
	return p3;
}

Term * Copy_Term (Term *t)
{
	Term *t3 = g_new0(Term, 1);
	GList *pn;
	for (pn = t->parts; pn; pn=pn->next)
	{
		Part *p = pn->data;
		t3->parts = g_list_append (t3->parts, Copy_Part(p));
	}

	t3->coeff = t->coeff;
	return t3;
}

Expr *Copy_Expr (Expr *e)
{
	Expr *e3 = g_new0(Expr, 1);

	GList *en;
	for (en = e->terms; en; en=en->next)
	{
		Term *tr = en->data;
		e3->terms = g_list_append (e3->terms, Copy_Term(tr));
	}
	return e3;
}

/* ======================================================== */

Expr *Add_E_T (Expr *e, Term *t)
{
	Expr *e3 = Copy_Expr(e);

	e3->terms = g_list_append (e3->terms, Copy_Term(t));
	return e3;
}

/* ======================================================== */
/* Derivatives */

Expr *DD_Term (Term *t, int which)
{
	if (0 == t->coeff) return NULL;
	if (NULL == t->parts) return NULL;

	Expr *e3 = g_new0(Expr, 1);

	GList *pn;
	for (pn = t->parts; pn; pn=pn->next)
	{
		Term *t3 = g_new0 (Term,1);
		t3->coeff = t->coeff;

		GList *pm;
		for (pm = t->parts; pm; pm=pm->next)
		{
			Part *p = pm->data;
			Part *p3 = Copy_Part (p);
			
			if (pm == pn)
			{
				if (1 == which) p3->dx ++;
				if (2 == which) p3->dy ++;
			}
			t3->parts = g_list_append (t3->parts, p3);
		}
		e3 = Add_E_T (e3, t3);		
	}

	return e3;
}

Expr *DDX_Term (Term *t)
{
	return DD_Term (t,1);
}

Expr *DDY_Term (Term *t)
{
	return DD_Term (t,2);
}

Expr *DD (Expr *e, int which)
{
	Expr *e3 = g_new0(Expr, 1);
	
	GList *en;
	for (en = e->terms; en; en=en->next)
	{
		Term *tr = en->data;
		Expr *em = DD_Term (tr, which);
		if (!em) continue;
		e3->terms = g_list_concat (e3->terms, em->terms);
		g_free (em);
	}
	return e3;
}

Expr *DDX (Expr *e)
{
	return DD (e,1);
}

Expr *DDY (Expr *e)
{
	return DD (e,2);
}


/* ======================================================== */
/* userland */

Term *tger (int x)
{
	Term *t = g_new (Term,1);
	t->coeff = x;
	return t;
}

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

/* ======================================================== */

/* return metrix xx */
Expr * gxx (void)
{
	Expr *e = g_new0 (Expr,1);
	e = Add_E_T (e, tger(1));
	e = Add_E_T (e, Mult_T_T (dxf(), dxf()));
	return e;
}

/* return metric yy */
Expr * gyy (void)
{
	Expr *e = g_new0 (Expr,1);
	e = Add_E_T (e, tger(1));
	e = Add_E_T (e, Mult_T_T (dyf(), dyf()));
	return e;
}

/* return metric xy */
Expr * gxy (void)
{
	Expr *e = g_new0 (Expr,1);
	e = Add_E_T (e, Mult_T_T (dxf(), dyf()));
	return e;
}

void setup (void)
{
	Expr *e = gxx();

	printf ("gxx is:\n");
	PrintExpr (e);

	e = gxy();
	printf ("gxy is:\n");
	PrintExpr (e);

	e = gyy();
	printf ("gyy is:\n");
	PrintExpr (e);

	printf ("dx gxx is:\n");
	PrintExpr (DDX(gxx()));
}

main ()
{

	setup();
}
