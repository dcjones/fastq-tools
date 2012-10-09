/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * parse :
 * A parser for FASTQ files.
 *
 */

#ifndef FASTQ_TOOLS_PARSE_H
#define FASTQ_TOOLS_PARSE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* A string structure to keep-track of a reserved space. */
typedef struct
{
    char*  s;    /* null-terminated string */
    size_t n;    /* length of s */
    size_t size; /* bytes allocated for s */
} str_t;


/* A single fastq entry. */
typedef struct
{
    str_t id1;
    str_t seq;
    str_t id2;
    str_t qual;
} seq_t;


/* Allocate a new empty seq_t. */
seq_t* seq_create();


/* Free a seq allocated with seq_create. */
void seq_free(seq_t* seq);


/* Hash a fastq entry. */
uint32_t seq_hash(const seq_t* seq);


/* Internal data for the fastq parser. */
typedef struct fastq_t_ fastq_t;


/* Create a new fastq parser object.
 *
 * Args:
 *   file: A file that has been opened for reading.
 */
fastq_t* fastq_create(FILE* file);


/* Free memory associated with a fastq_t object. */
void fastq_free(fastq_t*);


/* Read one fastq entry.
 *
 * Args:
 *   f: A fastq_t parser object.
 *   seq: A seq_t object that has been allocated with seq_create.
 *
 * Returns:
 *   True if an entry was read, false if end-of-file was reached.
 */
bool fastq_read(fastq_t* f, seq_t* seq);


/* Rewind the fastq file.
 *
 * The FILE passed to fastq_create must be seekable for this to work.
 */
void fastq_rewind(fastq_t* f);


/* Print a fastq entry. */
void fastq_print(FILE* fout, const seq_t* seq);


#endif

