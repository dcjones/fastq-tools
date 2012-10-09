
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-match :
 * Smith-Waterman alignments against sequences within a fastq file.
 *
 */


#include "common.h"
#include "parse.h"
#include "sw.h"
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif


static const char* prog_name = "fastq-match";


void print_help()
{
    fprintf(stdout, 
"fastq-match [OPTION]... QUERY [FILE]...\n"
"Perform Smith-Waterman local alignment of a query sequence\n"
"against each sequence in a fastq file.\n\n"
"Options:\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}


void fastq_match(FILE* fin, FILE* fout, sw_t* sw)
{
    int score;

    fastq_t* fqf = fastq_create(fin);
    seq_t* seq = seq_create();

    while (fastq_read(fqf, seq)) {
        fprintf(fout, "%s\t", seq->seq.s);

        fastq_sw_conv_seq((unsigned char*)seq->seq.s, seq->seq.n);
        score = fastq_sw(sw, (unsigned char*)seq->seq.s, seq->seq.n);

        fprintf(fout, "%d\n", score);
    }

    seq_free(seq);
    fastq_free(fqf);
}



int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    unsigned char* query;
    int query_len;

    sw_t* sw;

    FILE*  fin;

    int opt;
    int opt_idx;

    static struct option long_options[] =
        { 
          {"help",       no_argument,       NULL, 'h'},
          {"version",    no_argument,       NULL, 'V'},
          {0, 0, 0, 0}
        };


    while (1) {
        opt = getopt_long(argc, argv, "hV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                }
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


    if (optind >= argc) {
        fprintf(stderr, "A query sequence must be specified.\n");
        return 1;
    }

    query = (unsigned char*)argv[optind++];
    query_len = strlen((char*)query);
    fastq_sw_conv_seq(query, query_len);

    sw = fastq_alloc_sw(query, query_len);

    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        fastq_match(stdin, stdout, sw);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_match(fin, stdout, sw);
        }
    }

    fastq_free_sw(sw);

    return 0;
}



