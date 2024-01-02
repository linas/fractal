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
Open(const char *name, const char *ext)
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

FILE *Fopen(const char *name, const char *ext)
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

FILE *Fopenr(const char *name, const char *ext)
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

/* ----------------------- END OF FILE ---------------------- */
