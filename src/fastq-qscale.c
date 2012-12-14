/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-qualscale:
 * Determine the quality score encoding used by a fastq file.
 *
 */

#include <getopt.h>

#include "common.h"
#include "parse.h"

/* The scales here are taken from the wikipedia article for FASTQ, the relavent
 * part is reproduced here:
 *
 *  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 *  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
 *  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 *  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
 *  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
 *  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 *  |                         |    |        |                              |                     |
 * 33                        59   64       73                            104                   126
 *
 *    S - Sanger        Phred+33,  raw reads typically (0, 40)
 *    X - Solexa        Solexa+64, raw reads typically (-5, 40)
 *    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 *    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
 *       with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
 *       (Note: See discussion above).
 *        L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 */


typedef struct qual_scale_t_
{
    const char* description;
    char min_qual, max_qual;
} qual_scale_t;


/* When the scale is ambiguous, we choose the first compatible one. Hence these
 * are ordered roughly by increasing exclusivity. */
#define NUM_SCALES 5
static qual_scale_t scales[NUM_SCALES] =
{
    {"Sanger/Phred+33",       '!', 'I'},
    {"Illumina 1.8/Phred+33", '!', 'J'},
    {"Illumina 1.5/Phred+64", 'B', 'h'},
    {"Illumina 1.3/Phred+64", '@', 'h'},
    {"Solexa/Solexa+64",      ';', 'h'}
};


/* Return true if x has excatly one 1 bit. */
bool single_bit(uint32_t x)
{
    return x && !(x & (x - 1));
}


/* Make a bitset of compatible scales. */
uint32_t make_bitset(char min_qual, char max_qual)
{
    uint32_t s = 0;
    uint32_t i;
    for (i = 0; i < NUM_SCALES; ++i) {
        if (scales[i].min_qual <= min_qual &&
            scales[i].max_qual >= max_qual) {
            s |= 1 << i;
        }
    }

    return s;
}


void fastq_qualscale(const char* fn, FILE* fin)
{
    char min_qual = '~', max_qual = '!';

    fastq_t* fqf = fastq_create(fin);
    seq_t* seq = seq_create();

    /* Scales compatible with the data so far. */
    uint32_t compat_scales;

    size_t n = 100000;

    while (n-- && fastq_read(fqf, seq) && !single_bit(compat_scales)) {
        char* q = seq->qual.s;
        while (*q) {
            if (*q < min_qual) min_qual = *q;
            if (*q > max_qual) max_qual = *q;
            ++q;
        }

        compat_scales = make_bitset(min_qual, max_qual);
        if (compat_scales == 0 || single_bit(compat_scales)) break;
    }

    seq_free(seq);
    fastq_free(fqf);

    if (compat_scales == 0) {
        printf("%s: Unknown scale ['%c', '%c']\n", fn, min_qual, max_qual);
    }
    else {
        /* low order bit */
        unsigned int i;
        for (i = 0; !(compat_scales & (1 << i)); ++i) {}
        printf("%s: %s\n", fn, scales[i].description);
    }
}


static const char* prog_name = "fastq-qscale";


void print_help()
{
    fprintf(stderr,
"fastq-qscale [OPTION] [FILE]...\n"
"Detect and output the quality score scale for each file given as an argument.\n\n"
"Options:\n"
"  -h, --help              print this message\n"
"  -V, --version           output version information and exit\n"
    );
}


int main(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",    no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    int opt, opt_idx;
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

    if (optind >= argc) {
        fastq_qualscale("stdin", stdin);
    }
    else {
        for (; optind < argc; optind++) {
            FILE* fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_qualscale(argv[optind], fin);
            fclose(fin);
        }
    }

    return EXIT_SUCCESS;
}

