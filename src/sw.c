/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * This is a very simple, 'fast enough' implementation of the Smith-Waterman
 * algorithm specifically for short nucleotide sequences, working in O(mn) time
 * and O(m) space, implemented according to the original Gotoh paper and
 * Phil Green's implementation in cross_match.
 *
 * There is no fancy bit packing or vectorization, but such features would offer
 * diminishing returns when aligning short sequences such as high throughput
 * sequencing data. For example the Farrar SSE2 algorithm might be a tiny bit
 * faster, but would diminish portability.
 *
 */


#include "sw.h"
#include "common.h"
#include <string.h>
#include <stdlib.h>


static const int sw_default_d[25] =
    /* A   C   G   T   N */
    {  1, -2, -2, -2 , 0,
      -2,  1, -2, -2,  0,
      -2, -2,  1, -2,  0,
      -2, -2, -2,  1,  0,
       0,  0,  0,  0,  0  };

static const int sw_mat50_d[25] =
    {  2, -2,  0, -2 , 0,
      -2,  2, -2,  0,  0,
       0, -2,  2, -2,  0,
      -2,  0, -2,  2,  0,
       0,  0,  0,  0,  0  };

static const int sw_mat70_d[25] =
    {  2, -2, -1, -2 , 0,
      -2,  2, -2, -1,  0,
      -1, -2,  2, -2,  0,
      -2, -1, -2,  2,  0,
       0,  0,  0,  0,  0  };


static inline int imax(int a, int b)
{
    return a > b ? a : b;
}

static inline int imax4(int a, int b, int c, int d)
{
    return imax(imax(a,b), imax(c,d));
}



void fastq_sw_conv_seq(unsigned char* seq, int n)
{
    while (*seq && n) {
        switch (*seq) {
            case 'A' :
            case 'a':
            case 'U':
            case 'u':
                *seq = 0;
                break;

            case 'C':
            case 'c':
                *seq = 1;
                break;

            case 'G':
            case 'g':
                *seq = 2;
                break;

            case 'T':
            case 't':
                *seq = 3;
                break;

            case 'N':
            case 'n':
            default:
                *seq = 4;
        }

        seq++;
        n--;
    }
}


sw_t* fastq_alloc_sw(const unsigned char* subject, int size)
{
    sw_t* sw = malloc_or_die(sizeof(sw_t));

    sw->subject = malloc_or_die(size);
    memcpy(sw->subject, subject, size);

    /* default cost matrix */
    memcpy(sw->d, sw_default_d, 25 * sizeof(int)); 

    /* default gap costs */
    sw->gap_open   = -4;
    sw->gap_extend = -3;

    sw->work0 = malloc_or_die(size * sizeof(int));
    sw->work1 = malloc_or_die(size * sizeof(int));
    sw->size   = size;

    return sw;
}


void fastq_free_sw(sw_t* sw)
{
    free(sw->subject);
    free(sw->work0);
    free(sw->work1);
    free(sw);
}



int fastq_sw(sw_t* sw, const unsigned char* x, int n)
{
    /* conveniance */
    int*      maxstu = sw->work0;
    int*           t = sw->work1;
    int            m = sw->size;
    const int*     d = sw->d;
    int       gap_op = sw->gap_open;
    int       gap_ex = sw->gap_extend;
    unsigned char* y = sw->subject;

    int i, j;

    int score = 0;

    /* zero maxstu row */
    memset(maxstu, 0, m * sizeof(int));

    /* initialize t row to a negative value to prohibit
     * leading with gap extensions */
    for (j = 0; j < m; j++) t[j] = -1;

    int u, s;
    int maxstu0;

    for (i = 0; i < n; i++) {

        /* special case for column 0 */
        t[0] = imax(maxstu[0] + gap_op, t[0] + gap_ex);
        u    = gap_op;
        maxstu[0] = imax4(0, d[5 * y[0] + x[0]], t[0], u);
        maxstu0 = 0;
        score = imax(score, maxstu[0]);


        for (j = 1; j < m; j++) {
            t[j] = imax(maxstu[j] + gap_op, t[j] + gap_ex);
            u    = imax(maxstu[j-1] + gap_op, u + gap_ex);
            s    = maxstu0 + d[5 * y[j] + x[i]];
            maxstu0 = maxstu[j];

            maxstu[j] = imax4(0, s, t[j], u);
            score = imax(score, maxstu[j]);
        }
    }

    return score;
}



