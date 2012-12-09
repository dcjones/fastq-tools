/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * hash :
 * A quick and simple all-purpose hash table.
 *
 */


#ifndef FASTQ_TOOLS_HASH_H
#define FASTQ_TOOLS_HASH_H

#include <stdlib.h>
#include <stdint.h>


typedef struct hashed_value_
{
    char*    value;
    size_t   len;
    uint32_t count;
    struct hashed_value_* next;
} hashed_value;


typedef struct
{
    hashed_value** A; /* table proper */
    size_t n;         /* table size */
    size_t m;         /* hashed items */
    size_t max_m;     /* max hashed items before rehash */
} hash_table;


hash_table* create_hash_table();

void destroy_hash_table(hash_table*);

void inc_hash_table(hash_table*, const char* value, size_t len);

hashed_value** dump_hash_table(hash_table*);


#endif

