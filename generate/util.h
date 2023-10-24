/*
 * ancient code, 1987, 1989, 1993
 */


#include <stdio.h>

#ifdef  __cplusplus
extern "C" {
#endif

int Open(const char *name, const char *ext);
FILE *Fopen(const char *name, const char *ext);
FILE *Fopenr(const char *name, const char *ext);
void Size(int *width, int * height, const char * name, int fd, int bpp);
int Read(int fd, char *buf, int n);

#ifdef  __cplusplus
};
#endif
