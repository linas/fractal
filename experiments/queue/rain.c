/*
 * rainging down
 */

#define NMAX 1000000
signed char term[NMAX];
int path[NMAX/4];

void 
init_term (void)
{
	int i;
	for (i=0; i<NMAX; i++) term[i] = 0;
	term[0] = 1;
	term[1] = 1;
}

int 
is_term (int ns)
{
	int i, j;
	int ipath = 0;
	signed char result;

	// printf ("yo check=%d ======================================\n", ns);
	i = ns;
	while (1)
	{
		if (NMAX <= i) break;
		result = term[i];
		if ((0 < result) || (0 > result)) break;
		term[i] = -2;

		path[ipath] = i;
		ipath ++;
		if (NMAX/4 <= ipath) break;
		if (0 == i%2)
		{
			i = i/2;
		}
		else
		{
			// i = 3*i+1; descends very fast
			i = 3*i-1;
		}
		// printf ("yo ip=%d %d\n", ipath, i);
	}

	if (-2 == result) result = -1;
	if (NMAX/4 <= ipath) ipath = NMAX/4 -1;
	for (j=0; j<ipath; j++)
	{
		term[path[j]] = result;
	}
	if (NMAX/4 <= ipath) return 0;
	if (NMAX <= i) return 0;

	if (0 < result) return 1+ipath;
	if (0 > result) return -1;
	return 0;
}

int
main () 
{
	long nstart;
	long nterm;
	long nunknown;
	long ninf;
	
	ninf = 0;
	nterm = 0;
	nunknown = 0;
	init_term();
	printf ("done initing\n");
	for (nstart = 1; nstart < NMAX; nstart ++)
	{
		int rc;
		rc = is_term (nstart);
		
		if (0 < rc) nterm ++;
		else if (0 == rc) nunknown ++;
		else if (0 > rc) ninf ++;

		if (nstart %1000 == 0) 
		{
			printf ("term= %ld tot=%ld unk=%ld inf=%ld\n", nterm, nstart, nunknown, ninf);
		}
	}

}
