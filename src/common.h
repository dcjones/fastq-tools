/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * common :
 * A few common functions, primarily for crashing whilst retaining our dignity.
 *
 */

#ifndef FASTQ_TOOLS_COMMON_H
#define FASTQ_TOOLS_COMMON_H

#include <stdio.h>

void print_version(FILE* f, const char* prog_name);

void or_die(int b, const char* msg);

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char*, const char*);

/* Open a file for reading, creating it if it doesn't exist, and complaining if
 * it does. */
FILE* open_without_clobber(const char* filename);

#endif

