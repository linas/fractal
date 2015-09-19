/*
 * Calkin-Wilf tree fusc function
 *
 * September 2015
 */
#include <stdio.h>

unsigned long fusc(unsigned long n)
{
	if (0 == n) return 0;
	if (1 == n) return 1;
	unsigned long m = n >> 1;
	if (0 == n%2) return fusc(m);
	return fusc(m) + fusc(m+1);
}

int main(int argc, char * argv[])
{
	for (int i=0; i< 2048; i++)
	{
		printf("%lu	%lu\n", i, fusc(i));
	}

}
