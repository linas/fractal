
#include <stdio.h>
#include <stdlib.h>

int wincnt=0;
int failcnt=0;

int lam;

void check(int a, int b, int c, int d)
{
#define LAM lam
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
					// if (0 != smsa) continue;
					if (1 == smsb)
					{ found = true;

						double x = - ((double) d)/((double) c);
						double lo =  ((double) a)/((double) c);
						double hi =  ((double) b)/((double) d);
						printf("win (%d %d %d %d) (%d %d %d %d) (%d %d %d %d) %g = %g %g\n",
						a,b,c,d,p,q,r,s,smsa, smsb, smsc, smsd, x, lo, hi);
						int det = smsa*smsd - smsb*smsc;
						if (det != -1) {printf("aiiii %d\n", det); exit(1); }
						if (lo <= hi) { printf("wtf!!\n"); exit(1); }
						wincnt++;
						return;
					}
				}
			}
		}
	}
	if (!found)
	{
		double x = - ((double) d)/((double) c);
		double lo =  ((double) a)/((double) c);
		double hi =  ((double) b)/((double) d);
		printf("fail (%d %d %d %d) %g = %g %g \n",a,b,c,d, x, lo, hi);
		failcnt++;
		if (lo <= hi) { printf("wtf!!\n"); exit(1); }
		/* exit(1); */
	}
}


int
main(int argc, char * argv[])
{
	int lim = atoi(argv[1]);
	lam = atoi(argv[2]);
// #define LIM 6
#define LIM lim
	for(int a=-LIM; a<LIM; a++)
	{
		for (int b=-LIM; b<LIM; b++)
		{
			for (int c = -LIM; c<LIM; c++)
			{
				for (int d = -LIM; d<LIM; d++)
				{
					if (-1 != a*d-b*c) continue;
					double x = - ((double) d)/((double) c);
					if (x <= 0.0 || 1.0 < x) continue;
					check(a,b,c,d);
				}
			}
		}
	}
	double rat = ((double) wincnt)/ ((double) failcnt);
	printf("wincnct=%d %d %g\n", wincnt, failcnt, rat);
}
