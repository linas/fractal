/*
 * Copyright (c) 1987,1989 IBM Corporation
 * Copyright (c) 1989 Linas Vepstas
 * Author: Bruce Lucas
 * Author: Linas Vepstas
 */

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "util.h"

/*
 * Open with possible substitution of extension
 * name==0 => stdin
 */

int
Open(char *name, char *ext)
{
    int fd;

    if (ext) {
        /* substitute extension */
        char *slash, *dot, full[200];
        if (!name)
            return -1;
        strcpy(full, name);
        slash = (char *) rindex(full, '/');
        dot = (char *) rindex((slash?slash:full), '.');
        if (dot)
            *dot = '\0';
        strcat(full, ext);
        fd = open(full, O_RDONLY);
    } else if (!name) {
        /* stdin */
        return 0;
    } else {
        /* simple name */
        fd = open(name, O_RDONLY);
    }
    return fd;
}


/*
 * Linas Vepstas -- 6 July 1989
 * STDIO Fopen with possible substitution of extension
 * name==0 => stdin
 */

FILE *Fopen(char *name, char *ext)
{
    FILE *fyle;

    if (ext) {
        /* substitute extension */
        char *slash, *dot, full[200];
        if (!name)
            return NULL;
        strcpy(full, name);
        slash = (char *) rindex(full, '/');
        dot = (char *) rindex((slash?slash:full), '.');
        if (dot)
            *dot = '\0';
        strcat(full, ext);
        fyle = fopen (full, "wb");
    } else if (!name) {
        /* stdin */
        return NULL;
    } else {
        /* simple name */
        fyle = fopen (name, "wb");
    }
    return fyle;
}

/*
 * Linas Vepstas -- 6 July 1989
 * STDIO Fopen with possible substitution of extension
 * name==0 => stdin
 */

FILE *Fopenr(char *name, char *ext)
{
    FILE *fyle;

    if (ext) {
        /* substitute extension */
        char *slash, *dot, full[200];
        if (!name)
            return NULL;
        strcpy(full, name);
        slash = (char *) rindex(full, '/');
        dot = (char *) rindex((slash?slash:full), '.');
        if (dot)
            *dot = '\0';
        strcat(full, ext);
        fyle = fopen (full, "rb");
    } else if (!name) {
        /* stdin */
        return NULL;
    } else {
        /* simple name */
        fyle = fopen (name, "rb");
    }
    return fyle;
}


/*
 * Look up size in name.size,
 * or try to guess from the file size
 */

void
Size(int *width, int * height, char * name, int fd, int bpp)
{
    if (*width==0 || *height==0) {
        char buf[100];
        int fd = Open(name, ".size");
        int n = read(fd, buf, sizeof(buf));
        if (n >= 0) {
            buf[n] = 0;
            sscanf(buf, "%dx%d", width, height);
        }
    }
    if (*width==0 || *height==0) {
        int i, size = lseek(fd, 0, 2) / bpp;
        lseek(fd, 0, 0);
        for (i=1; i*i<=size; i++)
            if (size%i == 0)
                *height = i;
        if (*height) {
            *width = size / *height;
            fprintf(stderr, "guessing size is %dx%d\n", *width, *height);
        }
    }
    if (*width==0)
        fprintf(stderr, "please specify size");
        exit (1);
    if (*height==0)
        *height = *width;
}


/*
 * Try hard to read n bytes;
 * needed for input from pipes
 */

int
Read(int fd, char *buf, int n)
{
    int i, count;
    for (count=0; count<n; count+=i) {
        i = read(fd, buf+count, n-count);
        if (i < 0) {
            perror( "read");
            exit(1);
        }
        else if (i==0) {
            break;
        }
    }
    return count;
}

/* ----------------------- END OF FILE ---------------------- */
