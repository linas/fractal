
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

typedef struct eff_s Eff;

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


Expr * setup (void)
{
	Expr *e = g_new0 (Expr,0);

	return e;
}

void PrintExpr (Expr *e)
{
	printf ("Print expr: \n");
}

main ()
{

	Expr *e = setup();
}
