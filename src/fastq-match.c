
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-match :
 * Smith-Waterman alignments against sequences within a fastq file.
 *
 */


#include "fastq-common.h"
#include "fastq-parse.h"
#include "swsse2/blosum62.h"
#include "swsse2/swsse2.h"
#include "swsse2/matrix.h"
#include "swsse2/swstriped.h"
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


static int help_flag;


void print_help()
{
    fprintf( stderr, 
"fastq-match [OPTION]... QUERY [FILE]...\n"
"Perform Smith-Waterman local alignment of a query sequence\n"
"against each sequence in a fastq file.\n\n"
"Options:\n"
"  -h, --help              print this message\n"
    );
}


void convert_sequence(unsigned char* s, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        s[i] = (char)AMINO_ACID_VALUE[(int)s[i]];
    }
}


void fastq_match(FILE* fin, FILE* fout,
                 SwStripedData* sw_data,
                 unsigned char* query, int n,
                 SEARCH_OPTIONS* options)
{
    int score;

    int gap_init  = -(options->gapInit + options->gapExt);
    int gap_ext   = -options->gapExt;
    int threshold = options->threshold;

    fastq_t* fqf = fastq_open(fin);
    seq_t* seq = fastq_alloc_seq();

    while (fastq_next(fqf, seq)) {
        fprintf(fout, "%s\t", seq->seq.s);

        convert_sequence((unsigned char*)seq->seq.s, seq->seq.n);

        score = swStripedByte(query, n,
                              (unsigned char*)seq->seq.s, seq->seq.n,
                              gap_init, gap_ext,
                              sw_data->pvbQueryProf,
                              sw_data->pvH1,
                              sw_data->pvH2,
                              sw_data->pvE,
                              sw_data->bias);
        if (score >= 255) {
            score = swStripedWord(query, n,
                                  (unsigned char*)seq->seq.s, seq->seq.n,
                                  gap_init, gap_ext,
                                  sw_data->pvbQueryProf,
                                  sw_data->pvH1,
                                  sw_data->pvH2,
                                  sw_data->pvE);
        }

        fprintf(fout, "%d\n", score);
    }

    fastq_free_seq(seq);
    fastq_close(fqf);
}



int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    unsigned char* query;
    int query_len;
    SwStripedData* sw_data;
    signed char* mat = blosum62;
    SEARCH_OPTIONS options;

    options.gapInit   = -10;
    options.gapExt    = -2;
    options.threshold = -1;

    FILE*  fin;

    help_flag = 0;

    int opt;
    int opt_idx;

    static struct option long_options[] =
        { 
          {"help", no_argument, &help_flag, 1},
          {0, 0, 0, 0}
        };


    while (1) {
        opt = getopt_long(argc, argv, "h", long_options, &opt_idx);

        if( opt == -1 ) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                }
                break;

            case 'h':
                help_flag = 1;
                break;

            case '?':
                return 1;

            default:
                abort();
        }
    }

    if (help_flag) {
        print_help();
        return 0;
    }

    if (optind >= argc) {
        fprintf(stderr, "A query sequence must be specified.\n");
        return 1;
    }

    query = (unsigned char*)argv[optind++];
    query_len = strlen((char*)query);
    convert_sequence(query, query_len);

    sw_data = swStripedInit(query, query_len, mat);

    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        fastq_match(stdin, stdout, sw_data, query, query_len, &options);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_match(fin, stdout, sw_data, query, query_len, &options);
        }
    }

    swStripedComplete(sw_data);

    return 0;
}



