/*
 * ancient code, 1987, 1989, 1993
 */


#include <stdio.h>

#ifdef  __cplusplus
extern "C" {
#endif

int Open(char *name, char *ext);
FILE *Fopen(char *name, char *ext);
FILE *Fopenr(char *name, char *ext);
void Size(int *width, int * height, char * name, int fd, int bpp);
int Read(int fd, char *buf, int n);

#ifdef  __cplusplus
};
#endif
