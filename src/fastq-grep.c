

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


void print_help()
{
    fprintf( stderr, 
"fastq-grep [OPTION]... PATTERN [FILE]...\n"
"Search for PATTERN in the read sequence in each FILE or standard input.\n"
"PATTERN, by default, is a perl compatible regular expression."
    );

}

static int invert_flag;
static int help_flag;

int main(int argc, char* argv[])
{

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

    /* TODO */

    return 0;
}



