/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-sample :
 * Sample reads with or without replacement from a FASTQ file.
 *
 */

#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "common.h"
#include "parse.h"


typedef struct seq_array_t_
{
    seq_t* seqs;

    /* Number of seq objects. */
    size_t n;

    /* Size reserved in seqs. */
    size_t size;

    /* Space reserved for strings. */
    char* data;

    /* Data used. */
    size_t data_used;

    /* Total data size. */
    size_t data_size;
} seq_array_t;


seq_array_t* seq_array_create(size_t data_size)
{
    seq_array_t* a = malloc_or_die(sizeof(seq_array_t));
    a->size = 1024;
    a->n = 0;
    a->seqs = malloc_or_die(sizeof(seq_t));

    a->data_size = data_size;
    a->data_used = 0;
    a->data = malloc_or_die(data_size);

    return a;
}


void seq_array_free(seq_array_t* a)
{
    free(a->seqs);
    free(a->data);
    free(a);
}


/* Push a fastq entry to back of the array. Return false if there was not enough
 * space. */
bool seq_array_push(seq_array_t* a, const seq_t* seq)
{
    size_t size_needed = (seq->id1.n + 1) + (seq->seq.n + 1) +
                         (seq->id2.n + 1) + (seq->qual.n + 1);

    if (size_needed > a->data_size - a->data_used) return false;

    if (a->n == a->size) {
        a->size *= 2;
        a->seqs = realloc_or_die(a->seqs, a->size * sizeof(seq_t));
    }

    memcpy(seq->id1.s, &a->data[a->data_used], seq->id1.n + 1);
    a->seqs[a->n].id1.s = &a->data[a->data_used];
    a->seqs[a->n].id1.n = seq->id1.n + 1;
    a->data_used += seq->id1.n + 1;

    memcpy(seq->seq.s, &a->data[a->data_used], seq->seq.n + 1);
    a->seqs[a->n].seq.s = &a->data[a->data_used];
    a->seqs[a->n].seq.n = seq->seq.n + 1;
    a->data_used += seq->seq.n + 1;

    memcpy(seq->id2.s, &a->data[a->data_used], seq->id2.n + 1);
    a->seqs[a->n].id2.s = &a->data[a->data_used];
    a->seqs[a->n].id2.n = seq->id2.n + 1;
    a->data_used += seq->id2.n + 1;

    memcpy(seq->qual.s, &a->data[a->data_used], seq->qual.n + 1);
    a->seqs[a->n].qual.s = &a->data[a->data_used];
    a->seqs[a->n].qual.n = seq->qual.n + 1;
    a->data_used += seq->qual.n + 1;

    ++a->n;

    return true;
}


void seq_array_clear(seq_array_t* a)
{
    a->n = 0;
    a->data_used = 0;
}


void seq_array_sort(seq_array_t* a, int (*cmp)(const void*, const void*))
{
    qsort(a->seqs, a->n, sizeof(seq_t), cmp);
}


int seq_cmp_hash(const void* a_, const void* b_)
{
    const seq_t* a = (seq_t*) a_;
    const seq_t* b = (seq_t*) b_;
    /* TODO: hash and compare. */
    return 0;
}


static const char* prog_name = "fastq-sort";


void print_help()
{
    fprintf(stdout,
"fastq-sort [OPTION]... [FILE]...\n"
"Concatenate and sort FASTQ files and write to standard output.\n"
"Options:\n"
"  -h, --help      print this message\n"
"  -V, --version   output version information and exit\n"
   );
}



int main(int argc, char* argv[])
{
    int opt, opt_idx;
    size_t buffer_size = 100000000;
    int (*cmp)(const void*, const void*) = seq_cmp_hash;;

    static struct option long_options[] =
    {
        {"help",    no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    while (true) {
        opt = getopt_long(argc, argv, "hV", long_options, &opt_idx);
        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_help();
                return 0;

            case 'V':
                print_version(stdout, prog_name);
                return 0;

            case '?':
                return 1;

            default:
                abort();
        }
    }

    seq_array_t* a = seq_array_create(buffer_size);
    seq_t* seq = seq_create();

    fastq_t* f;
    if (optind >= argc) {
        f = fastq_create(stdin);
        while (fastq_read(f, seq)) {
            if (!seq_array_push(a, seq)) {
                seq_array_sort(a, cmp);

                /* TODO: dump a to a temporary file. Push that file name to an
                 * array somewhere.
                 * */

                seq_array_clear(a);
                if (seq_array_push(a, seq)) {
                    fprintf(stderr, "The buffer size is to small.\n");
                    return EXIT_FAILURE;
                }
            }
        }
        fastq_free(f);
    }
    else {
        FILE* file;
        for (; optind < argc; ++optind) {
            file = fopen(argv[optind], "rb");
            if (file == NULL) {
                fprintf(stderr, "Cannot open %s for reading.\n", argv[optind]);
                return EXIT_FAILURE;
            }
            f = fastq_create(file);

            while (fastq_read(f, seq)) {
                if (!seq_array_push(a, seq)) {
                    seq_array_sort(a, cmp);

                    /* TODO: dump a to a temporary file. Push that file name to
                     * an array somewhere.  */

                    seq_array_clear(a);
                    if (seq_array_push(a, seq)) {
                        fprintf(stderr, "The buffer size is to small.\n");
                        return EXIT_FAILURE;
                    }
                }
            }

            fastq_free(f);
            fclose(file);
        }
    }

    if (a->n > 0) {
        seq_array_sort(a, cmp);

        /* TODO: special case if everything fit into a. Just dump it to output.
         */

        /* TODO: dump to a temp file, push file name to the stack. */
    }

    /* TODO: Merge sort on the external files. */

    seq_free(seq);
    seq_array_free(a);

    return EXIT_SUCCESS;
}


