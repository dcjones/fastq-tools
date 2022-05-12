/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * fastq-grep :
 * Regular expression searches of the sequences within a FASTQ file.
 *
 */

#define _GNU_SOURCE

#include "common.h"
#include "parse.h"
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


static const char* prog_name = "fastq-tile-map";
static bool print_header = 0;

static float qsum[4][2][6][78] = {0};
static unsigned int nreads[4][2][6][78] = {0};
static unsigned int nperf[4][2][6][78] = {0};
static unsigned int nempty[4][2][6][78] = {0};
static unsigned int q20[4][2][6][78] = {0};
static unsigned int q30[4][2][6][78] = {0};

static float seq_mean_qual(const seq_t* s) {
  if (s->qual.n == 0) return 0.0;

  float q = 0.0;
  size_t i;
  for (i = 0; i < s->qual.n; ++i) {
    q += (float) s->qual.s[i];
  }
  return q / (float) s->qual.n;
}

bool allFs(char *s) {
    for (char *p = s; *p != 0; p++) {
      if (*p != 'F') {
        return false;
      }
    }

    return true;
}

bool allNs(char *s) {
    for (char *p = s; *p != 0; p++) {
      if (*p != 'N') {
        return false;
      }
    }

    return true;
}

void print_help() {
    fprintf(stdout,
"fastq-tile-map [FILE]...\n"
"Calculate several quality metrics for all image tiles referenced in any FILE or in standard input.\n"
"\n"
"Options:\n"
"  -H, --header        Output table header\n"
"  -h, --help          Print this message\n"
"  -V, --version       Print version information and exit\n"
    );
}

void scan_tiles(FILE* fin) {
  fastq_t* fqf = fastq_create(fin);
  seq_t* seq = seq_create();
  char colon[] = ":";

  while (fastq_read(fqf, seq)) {
    // seq->id1.s  # read name
    // seq->id1.n  # read name length
    // seq->seq.s  # sequence
    // seq->seq.n  # rsequenclength
    char *read_name = strdup(seq->id1.s);
    int lane = 0, surface = 0, swath = 0, tile = 0;  // allocation is important for sscanf()
    // fprintf(stderr, "read: %s\n", seq->id1.s);
    char *p = strtok(read_name, colon);
    int n = 1;  // Just scanned the first token
    while(p != NULL && ++n <= 5) {
      p = strtok(NULL, colon);
      // fprintf(stderr, "%d: %s\n", n, p);
      if (n == 4) {
        lane = *p - '0' - 1;
        // fprintf(stderr, "%d: lane %d\n", n, lane + 1);
      }
      if (n == 5) {
        sscanf(p, "%1d%1d%2d", &surface, &swath, &tile);
        --surface;
        --swath;
        --tile;

        // fprintf(stderr, "%d: %s -> %d, %d, %d\n", n, p, surface + 1, swath + 1, tile + 1);
      }
    }

    float q = seq_mean_qual(seq);
    qsum[lane][surface][swath][tile] += q;
    nreads[lane][surface][swath][tile] += 1;

    if ( q >= 20 + '!' ) {
      q20[lane][surface][swath][tile] += 1;
      if ( q >= 30 + '!' ) {
        q30[lane][surface][swath][tile] += 1;
      }
    }

    if (allFs(seq->qual.s)) {
      nperf[lane][surface][swath][tile] += 1;
    }
    if (allNs(seq->seq.s)) {
      nempty[lane][surface][swath][tile] += 1;
    }
  }

  seq_free(seq);
  fastq_free(fqf);
}



int main(int argc, char* argv[]) {
  SET_BINARY_MODE(stdin);


  int opt_idx;

  static struct option long_options[] = {
    {"help",         no_argument, NULL, 'h'},
    {"version",      no_argument, NULL, 'V'},
    {0, 0, 0, 0}
  };

  while (1) {
    int opt = getopt_long(argc, argv, "HhV", long_options, &opt_idx);

    if (opt == -1) break;

    switch (opt) {
      case 'h':
        print_help();
        return 0;

      case 'H':
        print_header = 1;
        break;

      case 'V':
        print_version(stdout, prog_name);
        return 0;

      case '?':
        return 1;

      default:
        abort();
    }

  }


  if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
    scan_tiles(stdin);
  }
  else {
    for (; optind < argc; optind++) {
      FILE*  fin = fopen(argv[optind], "rb");
      if (fin == NULL) {
        fprintf(stderr, "No such file '%s'.\n", argv[optind]);
        continue;
      }

      scan_tiles(fin);

      fclose(fin);
    }
  }

  if (print_header) {
    printf("Lane,Surface,Swath,Tile,Reads,Mean Q,Perfect Reads,Empty,Q20,Q30\n");
  }

  for (int lane = 0; lane < 4; lane++) {
    for (int surface = 0; surface < 2; surface++) {
      for (int swath = 0; swath < 6; swath++) {
        for (int tile = 0; tile < 78; tile++) {
          unsigned int n = nreads[lane][surface][swath][tile];
          if (n > 0) {
            printf(
              "%d,%d,%d,%d,%u,%5.2f,%u,%u,%u,%u\n",
              lane + 1,
              surface + 1,
              swath + 1,
              tile + 1,
              n,
              qsum[lane][surface][swath][tile] / n - '!',
              nperf[lane][surface][swath][tile],
              nempty[lane][surface][swath][tile],
              q20[lane][surface][swath][tile],
              q30[lane][surface][swath][tile]
            );
          }
          else {
            printf(
              "%d,%d,%d,%d,%u,%5.2f,%u,%u,%u,%u\n",
              lane + 1,
              surface + 1,
              swath + 1,
              tile + 1,
              (unsigned)0,
              0.0,
              (unsigned)0,
              (unsigned)0,
              (unsigned)0,
              (unsigned)0
            );
          }
        }
      }
    }
  }

  return 0;
}



