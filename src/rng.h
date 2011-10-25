/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * rng :
 * Robust pseudo-random number generation.
 *
 */

#ifndef FASTQ_TOOLS_RNG_H
#define FASTQ_TOOLS_RNG_H

typedef struct rng_t_ rng_t;

rng_t* fastq_rng_alloc();
void   fastq_rng_free(rng_t*);
void fastq_rng_seed(rng_t*, unsigned long seed);

/* Uniform integer in [0, k-1] */
unsigned long fastq_rng_uniform_int(rng_t*, unsigned long k);

#endif

