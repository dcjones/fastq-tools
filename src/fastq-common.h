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

void or_die(int b, const char* msg);

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char*, const char*);

#endif

