/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-parse :
 * A parser for FASTQ files.
 *
 * This parser is mostly derivative of Heng Li's.
 * See: http://lh3lh3.users.sourceforge.net/kseq.shtml
 *
 */

#ifndef FASTQ_TOOLS_PARSE_H
#define FASTQ_TOOLS_PARSE_H

#include <stdio.h>
#include <zlib.h>


typedef struct
{
    char*  s;    /* null-terminated string */
    size_t n;    /* length of s */
    size_t size; /* bytes allocated for s */
} str_t;



typedef struct
{
    str_t id1;
    str_t seq;
    str_t id2;
    str_t qual;
} seq_t;


seq_t* fastq_alloc_seq();
void fastq_free_seq(seq_t*);


typedef struct
{
    gzFile file;
    int    state;
    char*  buf;
    char*  c;
} fastq_t;


fastq_t* fastq_open(FILE*);
void fastq_close(fastq_t*);
int  fastq_next(fastq_t*, seq_t*);

#endif

