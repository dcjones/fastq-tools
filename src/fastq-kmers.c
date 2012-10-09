
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-kmers :
 * Tabulate k-mer frequencies with FASTQ files.
 *
 */

#include "common.h"
#include "parse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <zlib.h>


#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

static const char* prog_name = "fastq-kmers";

void print_help()
{
    fprintf( stderr, 
"fastq-kmers [OPTION]... [FILE]...\n"
"Print kmer counts for the given kmer size.\n"
"Output is in two tab-seperated columns for kmer and frequency.\n\n"
"Options:\n"
"  -k NUM, --size=NUM      kmer size (default: 1)\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}

static int k;

int packkmer( const char* s, uint32_t* kmer, int k )
{
    *kmer = 0;
    uint32_t nt;
    while (k--) {
        switch (s[k]) {
            case 'a':
            case 'A':
                nt = 0;
                break;
            case 'c':
            case 'C':
                nt = 1;
                break;
            case 'g':
            case 'G':
                nt = 2;
                break;
            case 't':
            case 'T':
                nt = 3;
                break;

            default:
                return 0;
        }

        *kmer = (*kmer << 2) | nt;
    }

    return 1;
}


void unpackkmer( uint32_t kmer, char* s, int k )
{
    int i;
    uint32_t nt;
    for (i = 0; i < k; i++, s++) {
        nt = kmer & 0x3;
        kmer = kmer >> 2;

        switch (nt) {
            case 0:
                *s = 'A';
                break;
            case 1:
                *s = 'C';
                break;
            case 2:
                *s = 'G';
                break;
            case 3:
                *s = 'T';
                break;
            default: break;
        }
    }

    *s = '\0';
}


void count_fastq_kmers(FILE* fin, uint32_t* cs)
{
    seq_t* seq = seq_create();
    fastq_t* fqf = fastq_create(fin);
    int i;
    int n;
    uint32_t kmer;

    while (fastq_read(fqf, seq)) {
        n = (int)seq->seq.n - k + 1;
        for (i = 0; i < n; i++) {
            if( packkmer(seq->seq.s + i, &kmer, k) ) {
                cs[kmer]++;
            }
        }
    }

    seq_free(seq);
    fastq_free(fqf);
}


void print_kmer_freqs(FILE* fout, uint32_t* cs)
{
    uint32_t n = 1 << (2*k); /* 4^k */

    char* kmer_str = malloc((k+1)*sizeof(char));
    uint32_t kmer;

    fprintf(fout, "kmer\tfrequency\n");
    for (kmer = 0; kmer < n; kmer++) {
        unpackkmer(kmer, kmer_str, k);
        fprintf(fout, "%s\t%u\n", kmer_str, cs[kmer]);
    }

    free(kmer_str);
}


int main(int argc, char* argv[])
{
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    k = 1;

    uint32_t n;   /* number of kmers: 4^k */
    uint32_t* cs; /* counts */

    FILE* fin;

    int opt;
    int opt_idx;
    static struct option long_options[] =
        { 
          {"size", no_argument,    0, 0},
          {"help", no_argument,    0, 'h'},
          {"version", no_argument, 0, 'V'},
          {0, 0, 0, 0}
        };

    while (1) {
        opt = getopt_long(argc, argv, "k:hV", long_options, &opt_idx);

        if( opt == -1 ) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                    if( opt_idx == 1 ) {
                        k = atoi(optarg);
                    }
                }
                break;

            case 'k':
                k = atoi(optarg);
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

    if (k < 1) {
        fprintf(stderr, "Kmer size must be at least 1.");
        return 1;
    }

    if (k > 16) {
        fprintf(stderr, "Kmer size must be at most 16.");
        return 1;
    }


    n = 1 << (2*k); /* i.e. 4^k */
    cs = malloc(n * sizeof(uint32_t));
    memset(cs, 0, n * sizeof(uint32_t));

    if (cs == NULL) {
        fprintf(stderr, "Insufficient memory to tally kmers of size %d\n", k );
        return 1;
    }

    if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
        count_fastq_kmers(stdin, cs);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            count_fastq_kmers(fin, cs);
        }
    }

    print_kmer_freqs( stdout, cs );
    free(cs);

    return 0;
}






