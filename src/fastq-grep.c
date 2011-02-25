

/* 
 * fastq-grep: regular expression searches of the sequences within a fastq file
 *
 * Febuary 2011 / Daniel Jones <dcjones@cs.washington.edu> 
 *
 */


#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>
#include <pcre.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


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
    );
}

static int invert_flag;
static int help_flag;
static int zout_flag;



void print_fastq_entry( FILE* fout, kseq_t* seq )
{
    fprintf(fout, "@%s\n%s\n+%s\n%s\n",
                  seq->name.s,
                  seq->seq.s,
                  seq->comment.s,
                  seq->qual.s );
}


void fastq_grep( gzFile fin, FILE* fout, pcre* re )
{
    int rc;
    int ovector[3];

    kseq_t* seq;
    seq = kseq_init(fin);

    while (kseq_read(seq) >= 0) {
        rc = pcre_exec(re,          /* pattern */
                       NULL,        /* extre data */
                       seq->seq.s,  /* subject */
                       seq->seq.l,  /* subject length */
                       0,           /* subject offset */
                       0,           /* options */
                       ovector,     /* output vector */
                       3         ); /* output vector length */

        if ((invert_flag && rc == PCRE_ERROR_NOMATCH) || rc >= 0) {
            print_fastq_entry( fout, seq );
        }
    }

    kseq_destroy(seq);
}



int main(int argc, char* argv[])
{
    const char* pat;
    pcre* re;
    const char* pat_error;
    int pat_error_offset;

    FILE*  fin;
    gzFile gzfin;


    invert_flag = 0;
    help_flag   = 0;
    zout_flag   = 0;

    int opt;
    int opt_idx;


    static struct option long_options[] =
        { 
          {"help", no_argument, &help_flag, 1},
          {"invert-match", no_argument, &invert_flag, 1},
          {0, 0, 0, 0}
        };

    while (1) {
        opt = getopt_long(argc, argv, "hv", long_options, &opt_idx);

        if( opt == -1 ) break;

        switch (opt) {
            case 0:
                if (long_options[opt_idx].flag != 0) break;
                if (optarg) {
                    /* TODO */
                }
                break;

            case 'h':
                help_flag = 1;
                break;

            case 'v':
                invert_flag = 1;
                break;

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
        gzfin = gzdopen( fileno(stdin), "rb" );
        if (gzfin == NULL) {
            fprintf(stderr, "Malformed file 'stdin'.\n");
            return 1;
        }

        fastq_grep(gzfin, stdout, re);

        gzclose(gzfin);
    }
    else {
        for (; optind < argc; optind++) {
            fin = fopen(argv[optind], "rb");
            if (fin == NULL) {
                fprintf(stderr, "No such file '%s'.\n", argv[optind]);
                continue;
            }

            gzfin = gzdopen(fileno(fin), "rb");
            if (gzfin == NULL) {
                fprintf(stderr, "Malformed file '%s'.\n", argv[optind]);
                continue;
            }

            fastq_grep(gzfin, stdout, re);

            gzclose(gzfin);
        }
    }

    pcre_free(re);

    return 0;
}



