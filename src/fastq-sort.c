/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-sort:
 * Sort fastq files efficiently.
 *
 */

#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

#include "common.h"
#include "parse.h"


/* User comparison function. */
static int (*user_cmp)(const void*, const void*);
static int (*cmp)(const void*, const void*);

int rev_cmp(const void* a, const void* b)
{
    return -user_cmp(a, b);
}


/* A collection of filenames of sorted chunks of fastq. */
typedef struct seq_dumps_t_
{
    char** fns;
    size_t n;
    size_t size;
} seq_dumps_t;


seq_dumps_t* seq_dumps_create()
{
    seq_dumps_t* d = malloc_or_die(sizeof(seq_dumps_t));
    d->n = 0;
    d->size = 64;
    d->fns = malloc_or_die(d->size * sizeof(char*));
    return d;
}


void seq_dumps_free(seq_dumps_t* d)
{
    size_t i;
    for (i = 0; i < d->n; ++i) {
        if (unlink(d->fns[i]) == -1) {
            fprintf(stderr, "Warning: unable to remove temporary file %s.\n",
                    d->fns[i]);
        }
        free(d->fns[i]);
    }
    free(d->fns);
    free(d);
}


static inline size_t heap_parent(size_t i) { return (i - 1) / 2; }
static inline size_t heap_left  (size_t i) { return 2 * i + 1; }
static inline size_t heap_right (size_t i) { return 2 * i + 2; }


/* Enqueue an item in the heap.
 *
 * Args:
 *   heap: A binary heap in which stores indexs into seqs.
 *   n: Fixed maximum size of the heap.
 *   m: Current size of the heap.
 *   seqs: Sequences.
 *   idx: Index to enqueue.
 *   cmp: Comparison function.
 *
 */
void heap_push(size_t* heap, size_t n, size_t* m, seq_t** seqs,
               int (*cmp)(const void*, const void*), size_t idx)
{
    if (*m >= n) return;
    size_t tmp, j, i = (*m)++;
    heap[i] = idx;

    /* percolate up */
    while (i > 0) {
        j = heap_parent(i);
        if (cmp(seqs[heap[i]], seqs[heap[j]]) < 0) {
            tmp = heap[i];
            heap[i] = heap[j];
            heap[j] = tmp;
            i = j;
        }
        else break;
    }
}


/* Dequeue an item from a heap. */
size_t heap_pop(size_t* heap, size_t* m, seq_t** seqs,
               int (*cmp)(const void*, const void*))
{
    assert(*m > 0);
    size_t ans = heap[0];
    heap[0] = heap[--(*m)];
    size_t tmp, l, r, j, i = 0;
    while (true) {
        l = heap_left(i);
        r = heap_right(i);

        if (l >= (*m)) break;

        if (r >= (*m)) {
            j = l;
        }
        else {
            j = cmp(seqs[heap[l]], seqs[heap[r]]) < 0 ? l : r;
        }

        if (cmp(seqs[heap[i]], seqs[heap[j]]) > 0) {
            tmp = heap[i];
            heap[i] = heap[j];
            heap[j] = tmp;
            i = j;
        }
        else break;
    }

    return ans;
}


/* n-way merge sort to stdout */
void merge_sort(const seq_dumps_t* d, int (*cmp)(const void*, const void*))
{
    FILE** files = malloc_or_die(d->n * sizeof(FILE*));
    size_t i;
    for (i = 0; i < d->n; ++i) {
        files[i] = fopen(d->fns[i], "rb");
        if (files[i] == NULL) {
            fprintf(stderr, "Cannot open temporary file %s for reading.\n",
                    d->fns[i]);
            exit(EXIT_FAILURE);
        }
    }

    fastq_t** fs = malloc_or_die(d->n * sizeof(fastq_t*));
    seq_t** seqs = malloc_or_die(d->n * sizeof(seq_t*));
    for (i = 0; i < d->n; ++i) {
        fs[i] = fastq_create(files[i]);
        seqs[i] = seq_create();
    }

    /* A binary heap of indexes to fs. We use this to repeatedly pop the
     * smallest fastq entry. */
    size_t* heap = malloc_or_die(d->n * sizeof(size_t));

    /* heap size */
    size_t m = 0;

    for (i = 0; i < d->n; ++i) {
        if (fastq_read(fs[i], seqs[i])) {
            heap_push(heap, d->n, &m, seqs, cmp, i);
        }
    }

    while (m > 0) {
        i = heap_pop(heap, &m, seqs, cmp);
        fastq_print(stdout, seqs[i]);
        if (fastq_read(fs[i], seqs[i])) {
            heap_push(heap, d->n, &m, seqs, cmp, i);
        }
    }

    for (i = 0; i < d->n; ++i) {
        seq_free(seqs[i]);
        fastq_free(fs[i]);
        fclose(files[i]);
    }

    free(files);
    free(fs);
}


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
    a->seqs = malloc_or_die(a->size * sizeof(seq_t));

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

    memcpy(&a->data[a->data_used], seq->id1.s, seq->id1.n + 1);
    a->seqs[a->n].id1.s = &a->data[a->data_used];
    a->seqs[a->n].id1.n = seq->id1.n + 1;
    a->data_used += seq->id1.n + 1;

    memcpy(&a->data[a->data_used], seq->seq.s, seq->seq.n + 1);
    a->seqs[a->n].seq.s = &a->data[a->data_used];
    a->seqs[a->n].seq.n = seq->seq.n + 1;
    a->data_used += seq->seq.n + 1;

    memcpy(&a->data[a->data_used], seq->id2.s, seq->id2.n + 1);
    a->seqs[a->n].id2.s = &a->data[a->data_used];
    a->seqs[a->n].id2.n = seq->id2.n + 1;
    a->data_used += seq->id2.n + 1;

    memcpy(&a->data[a->data_used], seq->qual.s, seq->qual.n + 1);
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


int seq_cmp_hash(const void* a, const void* b)
{
    uint32_t ha = seq_hash((seq_t*) a);
    uint32_t hb = seq_hash((seq_t*) b);

    if      (ha < hb) return -1;
    else if (ha > hb) return  1;
    else              return  0;
    return 0;
}


int seq_cmp_id(const void* a, const void* b)
{
    return strcmp(((seq_t*) a)->id1.s, ((seq_t*) b)->id1.s);
}


int seq_cmp_seq(const void* a, const void* b)
{
    return strcmp(((seq_t*) a)->seq.s, ((seq_t*) b)->seq.s);
}


void seq_array_dump(seq_dumps_t* d, const seq_array_t* a)
{
    const char* template = "/tmp/fastq_sort.XXXXXXXX";
    char* fn = malloc_or_die(strlen(template) + 1);
    memcpy(fn, template, strlen(template) + 1);
    if (mktemp(fn) == NULL) {
        fprintf(stderr, "Unable to create a temporary file.\n");
        exit(EXIT_FAILURE);
    }

    FILE* f = fopen(fn, "wb");
    if (f == NULL) {
        fprintf(stderr, "Unable to open temporary file %s for writing.\n", fn);
        exit(EXIT_FAILURE);
    }

    size_t i;
    for (i = 0; i < a->n; ++i) {
        fastq_print(f, &a->seqs[i]);
    }

    if (d->n == d->size) {
        d->size *= 2;
        d->fns = realloc_or_die(d->fns, d->size * sizeof(char*));
    }
    d->fns[d->n++] = fn;

    fclose(f);
}


static const char* prog_name = "fastq-sort";


void print_help()
{
    fprintf(stdout,
"fastq-sort [OPTION]... [FILE]...\n"
"Concatenate and sort FASTQ files and write to standard output.\n"
"Options:\n"
"  -I, --id           sort alphabetically by read identifier\n"
"  -S, --seq          sort alphabetically by sequence\n"
"  -R, --random       randomly shuffle the sequences\n"
"  -G, --gc           sort by GC content\n" /* TODO */
"  -M, --median-qual  sort by median quality score\n" /* TODO */
"  -h, --help         print this message\n"
"  -V, --version      output version information and exit\n"
   );
}


/* Parse a size specification, which is just a number with a K, M, G suffix. */
size_t parse_size(const char* str)
{
    char* endptr;
    unsigned long size = strtoul(str, &endptr, 10);

    if      (toupper(*endptr) == 'K') size *= 1000;
    else if (toupper(*endptr) == 'M') size *= 1000000;
    else if (toupper(*endptr) == 'G') size *= 1000000000;

    return size;
}


int main(int argc, char* argv[])
{
    int opt, opt_idx;
    size_t buffer_size = 100000000;
    bool reverse_sort = false;
    user_cmp = seq_cmp_id;

    static struct option long_options[] =
    {
        {"buffer-size", required_argument, NULL, 'S'},
        {"reverse",     no_argument, NULL, 'r'},
        {"id",          no_argument, NULL, 'i'},
        {"seq",         no_argument, NULL, 's'},
        {"random",      no_argument, NULL, 'R'},
        {"help",        no_argument, NULL, 'h'},
        {"version",     no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    while (true) {
        opt = getopt_long(argc, argv, "S:risRhV", long_options, &opt_idx);
        if (opt == -1) break;

        switch (opt) {
            case 'S':
                buffer_size = parse_size(optarg);
                break;

            case 'r':
                reverse_sort = true;
                break;

            case 'i':
                user_cmp = seq_cmp_id;
                break;

            case 's':
                user_cmp = seq_cmp_seq;
                break;

            case 'R':
                user_cmp = seq_cmp_hash;
                break;

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

    cmp = reverse_sort ? rev_cmp : user_cmp;

    seq_array_t* a = seq_array_create(buffer_size);
    seq_dumps_t* d = seq_dumps_create();
    seq_t* seq = seq_create();

    fastq_t* f;
    if (optind >= argc) {
        f = fastq_create(stdin);
        while (fastq_read(f, seq)) {
            if (!seq_array_push(a, seq)) {
                seq_array_sort(a, cmp);
                seq_array_dump(d, a);
                seq_array_clear(a);
                if (!seq_array_push(a, seq)) {
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
                    seq_array_dump(d, a);
                    seq_array_clear(a);
                    if (!seq_array_push(a, seq)) {
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

        /* We were able to sort everything in memory. */
        if (d->n == 0) {
            size_t i;
            for (i = 0; i < a->n; ++i) {
                fastq_print(stdout, &a->seqs[i]);
            }
        }
        else {
            seq_array_dump(d, a);
            merge_sort(d, cmp);
        }
    }

    seq_dumps_free(d);
    seq_free(seq);
    seq_array_free(a);

    return EXIT_SUCCESS;
}


