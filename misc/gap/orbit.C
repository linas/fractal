
/*
 * orbit.C
 *
 * Orbits of vectors under the dyadic representation 
 * of the modular group.
 *
 * Linas Vepstas November 2004
 */

double gee[2][2];
double are[2][2];

void 
init (void)
{
	gee[0][0] = 0.5;
	gee[0][1] = 0.0;
	gee[1][0] = 0.0;
	gee[1][1] = 1.0;

	are[0][0] = 1.0;
	are[0][1] = 0.0;
	are[1][0] = 1.0;
	are[1][1] = -1.0;
}

void
recur (double mat[2][2], int lvl, int width)
{
	double tmp[2][2];
	int i;
	lvl --;
	if (0 == lvl) return;

	for (i=0; i<width; i++)
	{
		tmp[0][0] = gee[0][0] * mat[0][0] + gee[0][1] * mat[1][0];
		tmp[0][1] = gee[0][0] * mat[0][1] + gee[0][1] * mat[1][1];
		tmp[1][0] = gee[1][0] * mat[0][0] + gee[1][1] * mat[1][0];
		tmp[1][1] = gee[1][0] * mat[0][1] + gee[1][1] * mat[1][1];

		mat[0][0] = tmp[0][0];
		mat[0][1] = tmp[0][1];
		mat[1][0] = tmp[1][0];
		mat[1][1] = tmp[1][1];

		printf ("%8.6g	%8.6g\n", mat[0][1], mat[1][1]);

		tmp[0][0] = are[0][0] * mat[0][0] + are[0][1] * mat[1][0];
		tmp[0][1] = are[0][0] * mat[0][1] + are[0][1] * mat[1][1];
		tmp[1][0] = are[1][0] * mat[0][0] + are[1][1] * mat[1][0];
		tmp[1][1] = are[1][0] * mat[0][1] + are[1][1] * mat[1][1];

		printf ("%8.6g	%8.6g\n", tmp[0][1], tmp[1][1]);

		recur (tmp, lvl, width);
	}
}

int main (int argc, char * argv[])
{

	double mat[2][2];
	
	mat[0][0] = 1.0;
	mat[0][1] = 0.0;
	mat[1][0] = 0.0;
	mat[1][1] = 1.0;

	init ();

	recur (mat, 3, 3);
}
