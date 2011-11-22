
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-qualadj:
 * Adjust quality scores by a given offset.
 *
 */


#include "common.h"
#include "parse.h"
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


static const char* prog_name = "fastq-grep";

void print_help()
{
    fprintf(stdout, 
"fastq-qualadj [OPTION]... OFFSET [FILE]...\n"
"The given offset is added to each and every quality score, where\n"
"the offset may be negative.\n"
"Options:\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}

void fastq_qualadj(FILE* fin, FILE* fout, int offset)
{
    fastq_t* fqf = fastq_open(fin);
    seq_t* seq = fastq_alloc_seq();
    size_t i;
    int c;

    while (fastq_next(fqf, seq)) {
        for (i = 0; i < seq->qual.n; ++i) {
            c = (int) seq->qual.s[i] - offset;
            c = c < 0 ? 0 : (c > 126 ? 126: c);
            seq->qual.s[i] = (char) c;
        }

        fastq_print(fout, seq);
    }

    fastq_free_seq(seq);
    fastq_close(fqf);
}


int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    char offset = 0;

    static struct option long_options[] =
        { 
          {"help",         no_argument, NULL, 'h'},
          {"version",      no_argument, NULL, 'V'},
          {0, 0, 0, 0}
        };

    int opt;
    int opt_idx;

    char* tmp;
    size_t tmplen;

    while (1) {
        opt = getopt_long(argc, argv, "hV0:1::2::3::4::5::6::7::8::9::", long_options, &opt_idx);

        if (opt == -1) break;

        switch(opt) {

            /* this is a bit of a hack to prevent getopt from choking on
             * negative numbers. */
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                tmplen = 2;
                if (optarg) tmplen += strlen(optarg);
                tmp = malloc(tmplen + 1);

                if (optarg) snprintf(tmp, tmplen + 1, "-%c%s", (char) opt, optarg);
                else        snprintf(tmp, tmplen + 1, "-%c",   (char) opt);

                offset = atoi(tmp);
                free(tmp);
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

    if (offset == 0 && optind >= argc) {
        fprintf(stderr, "An offset must be specified.\n");
        return 1;
    }
    else if (offset == 0) offset = atoi(argv[optind++]);

    FILE* fin;

    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        fastq_qualadj(stdin, stdout, offset);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_qualadj(fin, stdout, offset);

            fclose(fin);
        }
    }

    return 0;
}


