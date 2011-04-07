
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "common.h"
#include "version.h"
#include <stdlib.h>


void print_version(FILE *f, const char* prog_name)
{
    fprintf(f, "%s (fastq-tools) %s\n",
            prog_name, FASTQ_TOOLS_VERSION);
}


void or_die(int b, const char* msg)
{
    if (b == 0) {
        fputs(msg, stderr);
        exit(1);
    }
}


void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(1);
    }
    return p;
}


void* realloc_or_die(void* ptr, size_t n)
{
    void* p = realloc(ptr, n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(1);
    }
    return p;
}


FILE* fopen_or_die(const char* path, const char* mode)
{
    FILE* f = fopen(path, mode);
    if (f == NULL) {
        fprintf(stderr, "Can not open file %s with mode %s.\n", path, mode);
        exit(1);
    }
    return f;
}



