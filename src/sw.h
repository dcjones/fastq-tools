/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * sw :
 * Local alignments of nucleotide sequences via Smith-Waterman.
 *
 */

#ifndef FASTQ_TOOLS_SW_H
#define FASTQ_TOOLS_SW_H


typedef struct
{
    unsigned char* subject;
    int size;

    int d[25];      /* cost matrix */
    int gap_open;   /* gap open    */
    int gap_extend; /* gap extend  */

    /* matrix rows, used internally */
    int* work0;
    int* work1;

} sw_t;

/* convert a n ASCII nucleotide sequence to one suitable for fastq_sw */
void fastq_sw_conv_seq(unsigned char*, int n);

sw_t* fastq_alloc_sw(const unsigned char *subject, int size);
void  fastq_free_sw(sw_t*);

int fastq_sw(sw_t*, const unsigned char* query, int size);

#endif

