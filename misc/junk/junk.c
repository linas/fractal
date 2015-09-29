
#include <stdio.h>
#include <stdlib.h>

int wincnt=0;
int failcnt=0;

void check(int a, int b, int c, int d)
{
#define LAM 40
#define false 0
#define true 1
	int found = false;
	for(int p=-LAM; p<LAM; p++)
	{
		for (int r=-LAM; r<LAM; r++)
		{
			for (int q=-LAM; q<LAM; q++)
			{
				for (int s=-LAM; s<LAM; s++)
				{
					// int s = c*p-a*r;
					// int q = -(b*r + d*p);
					if (1 != p*s - q*r) continue;

					int smsa = p*a*s-c*p*q+r*s*b-r*q*d;
					int smsb = a*s*q-c*q*q+s*s*b-q*s*d;
					int smsc = -p*a*r+p*p*c-r*r*b+r*p*d;
					int smsd = -a*q*r+c*p*q-r*b*s+p*d*s;

					if (1 != smsc) continue;
					if (0 != smsd) continue;
					// if (1 == smsb) 
					{ found = true;

						double x = - ((double) d)/((double) c);
						printf("win (%d %d %d %d) (%d %d %d %d) (%d %d %d %d) %g\n",
						a,b,c,d,p,q,r,s,smsa, smsb, smsc, smsd, x);
						int det = smsa*smsd - smsb*smsc;
						if (det != 1) {printf("aiiii %d\n", det); exit(1); }
						wincnt++;
						return;
					}
				}
			}
		}
	}
	if (!found) { 
		double x = - ((double) d)/((double) c);
		printf("fail (%d %d %d %d) %g \n",a,b,c,d, x); failcnt++; /* exit(1); */}
}


int
main()
{
#define LIM 6
	for(int a=-LIM; a<LIM; a++)
	{
		for (int b=-LIM; b<LIM; b++)
		{
			for (int c = -LIM; c<LIM; c++)
			{
				for (int d = -LIM; d<LIM; d++)
				{
					if (1 != a*d-b*c) continue;
					double x = - ((double) d)/((double) c);
					if (x <= 0.0 || 1.0 < x) continue;
					check(a,b,c,d);
				}
			}
		}
	}
	printf("wincnct=%d %d\n", wincnt, failcnt);
}
