
#include <stdio.h>
#include "FareyTree.h"

main()
{
	FareyIterator fi;

	int i, j=0;
	int lastlvl = 0;
	for (i=0; i<300; i++)
	{
		int n, d;
		int loc = fi.GetNextFarey (&n, &d);
		int lvl = fi.GetLevel ();

		if (lvl != lastlvl) { printf ("\n"); lastlvl = lvl; j=0; }
		printf (" duude its %d  loc = %d  lvl=%d %d  %d/%d\n", i, loc, lvl, j, n, d);
		j++;
	}

}
