/*
 * Copyright (c) 1987,1989 IBM Corporation
 * Author: Bruce Lucas
 */

#include <fcntl.h>
#include <stdio.h>


/*
 * Open with possible substitution of extension
 * name==0 => stdin
 */

Open(name, ext)
char *name, *ext;
{
    int fd, i, dot, slash;

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

FILE *Fopen(name, ext)
char *name, *ext;
{
    int i, dot, slash;
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

FILE *Fopenr(name, ext)
char *name, *ext;
{
    int i, dot, slash;
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

static orient = 'l';

SetOrient(c)
char c;
{
    orient = c;
}

Size(width, height, name, fd, bpp)
char *name;
int fd;
int *width, *height;
int bpp;
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
            if (orient=='p') {
                /* use portrait mode */
                int temp;
                temp = *width; *width = *height; *height = temp;
            }
            fprintf(stderr, "guessing size is %dx%d\n", *width, *height);
        }
    }
    if (*width==0)
        Exit("please specify size");
    if (*height==0)
        *height = *width;
}


/*
 * Try hard to read n bytes;
 * needed for input from pipes
 */

Read(fd, buf, n)
int fd, n;
char *buf;
{
    int i, count;
    for (count=0; count<n; count+=i) {
        i = read(fd, buf+count, n-count);
        if (i < 0)
            Error("read");
        else if (i==0)
            break;
    }
    return count;
}


/*
 * Utilities
 */

Malloc(n)
{
    int m = malloc(n);
    if (!m)
        Exit("can't malloc");
    return m;
}

Exit(s)
char *s;
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}

Error(s)
char *s;
{
    perror(s);
    exit(1);
}




