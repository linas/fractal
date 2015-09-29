
#include <stdio.h>
#include <stdlib.h>

#define LIM 10
void check(int a, int b, int c, int d)
{
	for(int p=-LIM; p<LIM; p++)
	{
		for (int r=-LIM; r<LIM; r++)
		{
			int s = c*p-a*r;
			int q = -(b*r + d*p);

			if (1 != p*s - q*r) continue;

			int smsa = p*a*s-c*p*q+r*s*b-r*q*d;
			int smsb = a*s*q-c*q*q+s*s*b-q*s*d;
			int smsc = -p*a*r+p*p*c-r*r*b+r*p*d;
			int smsd = -q*q*r+c*p*q-r*b*s+p*d*s;

			printf("its (%d %d %d %d) (%d %d %d %d) (%d %d %d %d)\n",
			a,b,c,d,p,q,r,s,smsa, smsb, smsc, smsd);
			int det = smsa*smsd - smsb*smsc;
			if (det != 1) {printf("aiiii %d\n", det); exit(1); }
		}
	}
}


int
main()
{
	for(int a=-LIM; a<LIM; a++)
	{
		for (int b=-LIM; b<LIM; b++)
		{
			for (int c = -LIM; c<LIM; c++)
			{
				for (int d = -LIM; d<LIM; d++)
				{
					if (1 != a*d-b*c) continue;
					check(a,b,c,d);
				}
			}
		}
	}
}
