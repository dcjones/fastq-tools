
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-qual :
 * Collect quality score statistics.
 *
 */

#include "common.h"
#include "parse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>


#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

static const char* prog_name = "fastq-grep";

void print_help()
{
    fprintf(stdout, 
"fastq-qual [OPTION]... [FILE]...\n"
"Output a tab-delimnated table such that row i, column j, given \n"
"the number of times that quality score i occured in read position j\n\n."
"Options:\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n");
}



void tally_quals(FILE* fin, unsigned int** xs, size_t* n)
{
    seq_t* seq = fastq_alloc_seq();
    fastq_t* fqf = fastq_open(fin);

    size_t i;

    while (fastq_next(fqf, seq)) {
        if (seq->qual.n > *n) {
            *xs = realloc_or_die(*xs, 255 * seq->qual.n * sizeof(unsigned int));
            memset(*xs + *n, 0, 255 * (seq->qual.n - *n) * sizeof(unsigned int));
            *n  = seq->qual.n;
        }


        for (i = 0; i < seq->qual.n; ++i) {
            (*xs)[i * 255 + (int) seq->qual.s[i]] += 1;
        }
    }

    fastq_free_seq(seq);
    fastq_close(fqf);
}




void print_table(FILE* fout, unsigned int* xs, size_t n)
{
    size_t i, j;

    for (j = 0; j < 255; ++j) {
        fprintf(fout, "%u", xs[j]);
        for (i = 1; i < n; ++i) {
            fprintf(fout, "\t%u", xs[i * 255 + j]);
        }
        fputc('\n', fout);
    }
}



int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    FILE* fin;

    int opt;
    int opt_idx;
    static struct option long_options[] =
        { 
          {"help",    no_argument,    0, 'h'},
          {"version", no_argument, 0, 'V'},
          {0, 0, 0, 0}
        };

    while (1) {
        opt = getopt_long(argc, argv, "hV", long_options, &opt_idx);

        if( opt == -1 ) break;

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


    size_t n = 0;
    unsigned int* xs = NULL;

    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        tally_quals(stdin, &xs, &n);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            tally_quals(fin, &xs, &n);
        }
    }

    print_table(stdout, xs, n);

    free(xs);

    return 0;
}


