
#include <stdio.h>
#include "FareyTree.h"

main()
{
	FareyIterator fi;

	int i;
	for (i=0; i<300; i++)
	{
		int n, d;
		fi.GetNextFarey (&n, &d);

		printf (" duude its %d  %d/%d\n", i, n, d);
	}

}
