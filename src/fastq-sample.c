
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-sample :
 * Sample reads with or without replacement from a FASTQ file.
 *
 */

#include "common.h"
#include "parse.h"
#include "rng.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

static const char* prog_name = "fastq-sample";

static int replacement_flag;


void print_help()
{
    fprintf(stdout,
"fastq-sample [OPTION]... FILE [FILE2]\n"
"Sample random reads from a FASTQ file."
"Options:\n"
"  -n N                    the number of reads to sample (default: 10000)\n"
"  -p N                    the proportion of the total reads to sample\n"
"  -o, --output=PREFIX     output file prefix\n"
"  -c, --complement-output=PREFIX\n"
"                          output reads not included in the random sample to\n"
"                          a file (or files) with the given prefix (by default,\n"
"                          they are not output).\n"
"  -r, --with-replacement  sample with relpacement\n"
"  -s, --seed=SEED         a manual seed to the random number generator\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}


/* count the number of entries in a fastq file */
unsigned long count_entries(fastq_t* fqf)
{
    seq_t* seq = seq_create();
    unsigned long n = 0;
    while (fastq_read(fqf, seq)) ++n;
    seq_free(seq);

    return n;
}


/* compare two unsigned integers (for qsort) */
int cmpul( const void* a, const void* b ) {
    if(      *(unsigned long*) a < *(unsigned long*) b ) return -1;
    else if( *(unsigned long*) a > *(unsigned long*) b ) return 1;
    else                                                 return 0;
}


/* randomly shuffle an array of unsigned longs */
void shuffle(rng_t* rng, unsigned long* xs, unsigned long n)
{
    unsigned long i, j, k;
    for (i = n - 1; i > 0; --i) {
        j = fastq_rng_uniform_int(rng, i + 1);
        k = xs[j]; xs[j] = xs[i]; xs[i] = k; // swap
    }
}


unsigned long* index_with_replacement(rng_t* rng, unsigned long n, unsigned long k)
{
    unsigned long* xs = malloc_or_die(k * sizeof(unsigned long));
    size_t i;
    for (i = 0; i < k; ++i) xs[i] = fastq_rng_uniform_int(rng, n);
    return xs;
}


unsigned long* index_without_replacement(rng_t* rng, unsigned long n)
{
    unsigned long* xs = malloc_or_die(n * sizeof(unsigned long));
    size_t i;
    for (i = 0; i < n; ++i) xs[i] = i;
    shuffle(rng, xs, n);
    return xs;
}


void fastq_sample(unsigned long rng_seed,
                  const char* prefix, const char* cprefix,
                  FILE* file1, FILE* file2, unsigned long k, double p)
{
     /*
      * The basic idea is this:
      *
      * 1. Count the number of lines in the file, n.
      *
      * 2a. If sampling with replacement, generate k random integers in [0, n-1].
      *
      * 2b. If sampling without replacement, generate a list of integers 0..(n-1),
      *     shuffle with fisher-yates, then consider the first k.
      *
      * 3. Sort the integer list.
      *
      * 3. Read through the file again, when the number at the front of the integer
      *    list matches the index of the fastq etry, print the entry, and pop the
      *    number.
      */


    unsigned long n, n2;

    fastq_t* f1 = fastq_create(file1);
    fastq_t* f2 = file2 == NULL ? NULL : fastq_create(file2);

    n = count_entries(f1);
    if (f2 != NULL) {
        n2 = count_entries(f2);
        if (n != n2) {
            fprintf(stderr, "Input files have differing numbers of entries (%lu != %lu).\n", n, n2);
            exit(1);
        }
    }

    fastq_rewind(f1);
    if (f2 != NULL) fastq_rewind(f2);

    if (p > 0.0) {
        k = (unsigned long) (p * (double) n);
        if (!replacement_flag && k > n) k = n;
    }

    rng_t* rng = fastq_rng_alloc();
    fastq_rng_seed(rng, rng_seed);

    unsigned long* xs;
    if (replacement_flag) xs = index_with_replacement(rng, n, k);
    else                  xs = index_without_replacement(rng, n);

    qsort(xs, k, sizeof(unsigned long), cmpul);


    /* open output */
    FILE* fout1;
    FILE* fout2;

    char* output_name;
    size_t output_len;
    if (file2 == NULL) {
        output_len = strlen(prefix) + 7;
        output_name = malloc_or_die((output_len + 1) * sizeof(char));

        snprintf(output_name, output_len, "%s.fastq", prefix);
        fout1 = fopen(output_name, "wb");
        if (fout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        fout2 = NULL;

        free(output_name);
    }
    else {
        output_len = strlen(prefix) + 9;
        output_name = malloc_or_die((output_len + 1) * sizeof(char));

        snprintf(output_name, output_len, "%s.1.fastq", prefix);
        fout1 = fopen(output_name, "wb");
        if (fout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        snprintf(output_name, output_len, "%s.2.fastq", prefix);
        fout2 = fopen(output_name, "wb");
        if (fout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        free(output_name);
    }


    /* open complement output */
    FILE* cfout1 = NULL;
    FILE* cfout2 = NULL;

    if (cprefix != NULL && file2 == NULL) {
        output_len = strlen(cprefix) + 7;
        output_name = malloc_or_die((output_len + 1) * sizeof(char));

        snprintf(output_name, output_len, "%s.fastq", cprefix);
        cfout1 = fopen(output_name, "wb");
        if (cfout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        cfout2 = NULL;

        free(output_name);
    }
    else if (cprefix != NULL) {
        output_len = strlen(cprefix) + 9;
        output_name = malloc_or_die((output_len + 1) * sizeof(char));

        snprintf(output_name, output_len, "%s.1.fastq", cprefix);
        cfout1 = fopen(output_name, "wb");
        if (cfout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        snprintf(output_name, output_len, "%s.2.fastq", cprefix);
        cfout2 = fopen(output_name, "wb");
        if (cfout1 == NULL) {
            fprintf(stderr, "Cannot open file %s for writing.\n", output_name);
            exit(1);
        }

        free(output_name);
    }



    unsigned long i = 0; // read number
    unsigned long j = 0; // index into xs

    int ret;
    seq_t* seq1 = seq_create();
    seq_t* seq2 = seq_create();

    while (j < k && fastq_read(f1, seq1)) {
        if (f2 != NULL){
            ret = fastq_read(f2, seq2);
            if (ret == 0) {
                fputs("Input files have differing numbers of entries.\n", stderr);
                exit(1);
            }
        }

        if (xs[j] == i) {
            while (j < k && xs[j] == i) {
                fastq_print(fout1, seq1);
                if (f2 != NULL) fastq_print(fout2, seq2);
                ++j;
            }
        }
        else if (cfout1 != NULL) {
            fastq_print(cfout1, seq1);
            if (f2 != NULL) fastq_print(cfout2, seq2);
        }

        ++i;
    }

    seq_free(seq1);
    seq_free(seq2);
    fastq_free(f1);
    if (f2 != NULL) fastq_free(f2);

    fclose(fout1);
    if (fout2 != NULL) fclose(fout2);

    if (cfout1 != NULL) fclose(cfout1);
    if (cfout2 != NULL) fclose(cfout2);

    fastq_rng_free(rng);
    free(xs);
}


int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    int opt;
    int opt_idx;

    const char* prefix = NULL;
    const char* cprefix = NULL;
    unsigned long rng_seed = 4357;
    unsigned long k = 10000; // number of reads to sample
    double        p = -1;    // proportion of reads to sample


    static struct option long_options[] =
        { 
          {"with-replacement",  no_argument,       NULL, 'r'},
          {"complement-output", required_argument, NULL, 'c'},
          {"seed",              required_argument, NULL, 's'},
          {"output",            no_argument,       NULL, 'o'},
          {"help",              no_argument,       NULL, 'h'},
          {"version",           no_argument,       NULL, 'V'},
          {0, 0, 0, 0}
        };

    replacement_flag = 0;

    while (1) {
        opt = getopt_long(argc, argv, "n:p:o:rs:hV", long_options, &opt_idx);

        if( opt == -1 ) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                }
                break;

            case 'n':
                k = strtoul(optarg, NULL, 10);
                break;

            case 'p':
                p = atof(optarg);
                if (p < 0.0) {
                    fputs("Sample proportion ('-p') is less than zero.\n", stderr);
                    return 1;
                }
                break;

            case 'r':
                replacement_flag = 1;
                break;

            case 's':
                rng_seed = strtoul(optarg, NULL, 10);
                break;

            case 'o':
                prefix = optarg;
                break;

            case 'c':
                cprefix = optarg;
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

    FILE* file1 = NULL;
    FILE* file2 = NULL;

    char* prefix_alloc = NULL;

    if (optind >= argc) {
        fputs("An input file must be given.\n", stderr);
        print_help();
        exit(1);
    }
    else {
        file1 = fopen(argv[optind], "rb");
        if (file1 == NULL) {
            fprintf(stderr, "Cannot open '%s' for reading.\n", argv[optind]);
            return 1;
        }

        if (prefix == NULL) {
            /* guess at a reasonable output refix by trimming the
             * trailing file extension, if any.  */
            char* tmp;

            /* base name */
            tmp = strrchr(argv[optind], '/');
            if (tmp != NULL) argv[optind] = tmp + 1;

            /* exclude file suffixes */
            tmp = strchr(argv[optind], '.');
            if (tmp == NULL) prefix = argv[optind];
            else {
                prefix_alloc = malloc_or_die((tmp - argv[optind] + 1) * sizeof(char));
                memcpy(prefix_alloc, argv[optind], (tmp - argv[optind]) * sizeof(char));
                prefix_alloc[tmp - argv[optind]] = '\0';
                prefix = prefix_alloc;
            }
        }
        ++optind;

        if (optind < argc) {
            file2 = fopen(argv[optind], "rb");
            if (file2 == NULL) {
                fprintf(stderr, "Cannot open '%s' for reading.\n", argv[optind]);
                return 1;
            }
        }
    }

    fastq_sample(rng_seed, prefix, cprefix, file1, file2, k, p);

    free(prefix_alloc);
    return 0;
}








