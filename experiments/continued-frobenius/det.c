
typedef int vector[30];

void
alt (vector *ir, vector *ic, int vlen)
{
	if (2 == vlen)
	{
		int ra, rb, ca, cb;
		ra = (*ir)[0];
		rb = (*ir)[1];

		ca = (*ic)[0];
		cb = (*ic)[1];

		printf ("(((*m)[r%d][c%d]) * ((*m)[r%d][c%d]) - ((*m)[r%d][c%d]) * ((*m)[r%d][c%d]))\n",
			ra,ca,rb,cb, ra,cb,rb,ca);
	}
	else
	{
		vector row, col;

		int rrr = (*ir)[0];
		int i,j;
		for (i=0; i<vlen-1; i++)
		{
			row[i] = (*ir)[i+1];
		}

		int sign = 1;
		printf ("(  ");
		for (i=0; i<vlen; i++)
		{
			int k=0;
			for (j=0; j<vlen-1; j++)
			{
				if (k==i) k++;  // skip this column
				col[j] = (*ic)[k];
				k++;
			}

			int ccc = (*ic)[i];
			if (i != 0) {
				if (sign > 0) printf (" + ");
				else printf (" - ");
			}
			printf ("((*m)[r%d][c%d]) * ", rrr, ccc);
			alt (&row, &col, vlen-1);
			sign = -sign;
		}
		printf (")\n");
	}
}

void
determinant (int dim)
{
	int i;
	vector row, col;

	for (i=0; i<dim; i++)
	{
		row[i] = i;
		col[i] = i;
	}

	alt (&row, &col, dim);
}

main ()
{
	determinant (5);
	return 0;
}
