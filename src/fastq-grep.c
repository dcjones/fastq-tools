
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-grep :
 * Regular expression searches of the sequences within a FASTQ file.
 *
 */


#include "common.h"
#include "parse.h"
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>
#include <pcre.h>


#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif


void print_help()
{
    fprintf( stderr, 
"fastq-grep [OPTION]... PATTERN [FILE]...\n"
"Search for PATTERN in the read sequences in each FILE or standard input.\n"
"PATTERN, by default, is a perl compatible regular expression.\n\n"
"Options:\n"
"  -h, --help              print this message\n"
"  -v, --invert-match      select nonmatching entries\n"
"  -c, --count             output only the number of matching sequences\n"
    );
}

static int invert_flag;
static int help_flag;
static int count_flag;



void fastq_grep(FILE* fin, FILE* fout, pcre* re)
{
    int rc;
    int ovector[3];
    size_t count = 0;

    fastq_t* fqf = fastq_open(fin);
    seq_t* seq = fastq_alloc_seq();

    while (fastq_next(fqf, seq)) {
        rc = pcre_exec(re,          /* pattern */
                       NULL,        /* extre data */
                       seq->seq.s,  /* subject */
                       seq->seq.n,  /* subject length */
                       0,           /* subject offset */
                       0,           /* options */
                       ovector,     /* output vector */
                       3         ); /* output vector length */

        if ((invert_flag && rc == PCRE_ERROR_NOMATCH) || rc >= 0) {
            if (count_flag) count++;
            else            fastq_print(fout, seq);
        }
    }

    fastq_free_seq(seq);
    fastq_close(fqf);

    if (count_flag) fprintf(fout, "%zu\n", count);
}



int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    const char* pat;
    pcre* re;
    const char* pat_error;
    int pat_error_offset;

    FILE*  fin;


    invert_flag  = 0;
    help_flag    = 0;
    count_flag   = 0;

    int opt;
    int opt_idx;


    static struct option long_options[] =
        { 
          {"help", no_argument, &help_flag, 1},
          {"invert-match", no_argument, &invert_flag, 1},
          {"count", no_argument, &count_flag, 1},
          {0, 0, 0, 0}
        };

    while (1) {
        opt = getopt_long(argc, argv, "hvc", long_options, &opt_idx);

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

            case 'v':
                invert_flag = 1;
                break;

            case 'c':
                count_flag = 1;
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
        fprintf(stderr, "A pattern must be specified.\n");
        return 1;
    }

    pat = argv[optind++];
    re = pcre_compile( pat, PCRE_CASELESS, &pat_error, &pat_error_offset, NULL );


    if (re == NULL) {
        fprintf(stderr, "Syntax error in PCRE pattern at offset: %d: %s\n",
                pat_error_offset, pat_error );
        return 1;
    }


    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        fastq_grep(stdin, stdout, re);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_grep(fin, stdout, re);

            fclose(fin);
        }
    }

    pcre_free(re);

    return 0;
}



