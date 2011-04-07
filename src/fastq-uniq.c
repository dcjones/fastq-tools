
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-uniq :
 * Collapsing a fastq file into only unique read sequences.
 */


#include "common.h"
#include "hash.h"
#include "parse.h"
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

static const char* prog_name = "fastq-uniq";

void print_help()
{
    fprintf(stderr,
"fastq-uniq [OPTION] [FILE]...\n"
"Output a non-redundant FASTQ file, in which there are no duplicate reads.\n"
"(Warning: this program can be somewhat memory intensive.)\n\n"
"Options:\n"
"  -v, --verbose           print status along the way\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}


static int verbose_flag;
static size_t total_reads;



void fastq_hash(FILE* fin, hash_table* T)
{
    fastq_t* fqf = fastq_open(fin);
    seq_t* seq = fastq_alloc_seq();

    while (fastq_next(fqf, seq)) {
        inc_hash_table(T, seq->seq.s, seq->seq.n);

        total_reads++;
        if (verbose_flag && total_reads % 100000 == 0) {
            fprintf(stderr, "%zu reads processed ...\n", total_reads);
        }
    }

    fastq_free_seq(seq);
    fastq_close(fqf);
}


int compare_hashed_value_count(const void* x, const void* y)
{
    hashed_value* const * a = x;
    hashed_value* const * b = y;

    if( (*a)->count > (*b)->count ) return -1;
    if( (*a)->count < (*b)->count ) return 1;
    return 0;
}



void print_hash_table(FILE* fout, hash_table* T)
{
    hashed_value** S = dump_hash_table(T);
    qsort(S, T->m, sizeof(hashed_value*), compare_hashed_value_count);

    size_t i;
    for (i = 0; i < T->m; i++) {
        fprintf(fout, ">unique-read-%07zu (%zu copies)\n", i, S[i]->count);
        fwrite(S[i]->value, S[i]->len, sizeof(char), fout);
        fprintf(fout, "\n");
    }
    free(S);
}



int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    hash_table* T = create_hash_table();

    FILE* fin   ;

    int opt;
    int opt_idx;

    static struct option long_options[] =
    {
        {"verbose", no_argument, &verbose_flag, 1},
        {"help",    no_argument, NULL,          'h'},
        {"version", no_argument, NULL,          'V'},
        {0, 0, 0, 0}
    };

    while (1) {
        opt = getopt_long(argc, argv, "vhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                }
                break;

            case 'v':
                verbose_flag = 1;
                break;

            case '?':
                return 1;

            case 'h':
                print_help();
                return 0;

            case 'V':
                print_version(stdout, prog_name);
                return 0;

            default:
                abort();
        }
    }


    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        fastq_hash(stdin, T);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_hash(fin, T);
        }
    }

    print_hash_table(stdout, T);

    destroy_hash_table(T);
    return 0;
}
