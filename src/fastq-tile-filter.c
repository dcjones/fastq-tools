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


static const char* prog_name = "fastq-tile-filter";

static bool filter[4][2][6][78] = {0};

void print_help() {
    fprintf(stdout,
"\n"
"fastq-tile-filter [OPTIONS] [LIST] [FILE]\n"
"\n"
"Use the list of tiles to filter the input file or stdin.\n"
"\n"
"The comma-separated list must contain these columns in the first four positions:\n"
"\r"
"   Lane,Surface,Swath,Tile\n"
"\n"
"The header is optional.\n"
"\n"

"Options:\n"
"  -h, --help          Print this message\n"
"  -V, --version       Print version information and exit\n"
"  -v, --invert        Output reads associated with tiles that are not on the list\n"
    );
}

static int invert_flag;


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

        // fprintf(stderr, "%d: %s -> %d, %d, %d -> %d\n", n, p, surface + 1, swath + 1, tile + 1, filter[lane][surface][swath][tile]);
        if (filter[lane][surface][swath][tile]) {
          fastq_print(stdout, seq);
        }
      }
    }
  }

  seq_free(seq);
  fastq_free(fqf);
}



int main(int argc, char* argv[]) {
  SET_BINARY_MODE(stdin);
  SET_BINARY_MODE(stdout);
  const char* list_file;


  invert_flag = 0;
  int opt_idx;

  static struct option long_options[] = {
    {"help",         no_argument, NULL, 'h'},
    {"version",      no_argument, NULL, 'V'},
    {"invert",       no_argument, &invert_flag, 1},
    {0, 0, 0, 0}
  };

  while (1) {
    int opt = getopt_long(argc, argv, "HhVv", long_options, &opt_idx);

    if (opt == -1) break;

    switch (opt) {
      case 'h':
        print_help();
        return 0;

      case 'V':
        print_version(stdout, prog_name);
        return 0;

      case 'v':
        invert_flag = 1;
        for (int lane = 0; lane < 4; lane++) {
          for (int surface = 0; surface < 2; surface++) {
            for (int swath = 0; swath < 6; swath++) {
              for (int tile = 0; tile < 78; tile++) {
                filter[lane][surface][swath][tile] = 1;
              }
            }
          }
        }

        break;

      case '?':
        return 1;

      default:
        abort();
    }

  }

  if (optind >= argc) {
    fprintf(stderr, "The list of tiles to filter must be provided.\n");
    return 1;
  }

  list_file = argv[optind++];

  unsigned char magic[2];

  if (optind >= argc || (argc - optind == 1 && strcmp(argv[optind],"-") == 0)) {
    char buffer[1000];
    FILE*  fl = fopen(list_file, "rb");
    if (fl == NULL) {
      fprintf(stderr, "Could not open %s\n", argv[optind]);
      return 1;
    }

    // Parse the tile list
    int line = 0;
    while (fgets(buffer, 1000, fl)) {
      ++line;
      char *token = strtok(buffer, ",");
      int tn = 0;
      unsigned short lane = 0, surface = 0, swath = 0, tile = 0;

      while (token) {
        ++tn;
        // printf("line %d, token %d: %s\n", line, tn, token);

        // Skip the header, but check that it contains the necessary columns if present
        switch (tn) {
          case 1:
            if (!sscanf(token, "%hu", &lane)) {
              if (line == 1) { // Check if it is the lane header
                if (!(strncmp(token, "\"Lane\"", 6) == 0 || strncmp(token, "\"lane\"", 6) == 0 || strncmp(token, "Lane", 4) == 0 || strncmp(token, "lane", 4) == 0)) {
                  fprintf(stderr, "The first column in the filter list must be lane number; got %s instead", token);
                  return 1;
                }
              }
              else {
                fprintf(stderr, "%s:%d: Not an integer: %s", list_file, line, token);
                return 1;
              }
            }
            break;

          case 2:
            if (!sscanf(token, "%hu", &surface)) {
              if (line == 1) { // Check if it is the surface header
                if (!(strncmp(token, "\"Surface\"", 9) == 0 || strncmp(token, "\"surface\"", 9) == 0 || strncmp(token, "Surface", 7) == 0 || strncmp(token, "surface", 7) == 0)) {
                  fprintf(stderr, "The second column in the filter list must be surface number; got %s instead", token);
                  return 1;
                }
              }
              else {
                fprintf(stderr, "%s:%d: Not an integer: %s", list_file, line, token);
                return 1;
              }
            }
            break;

          case 3:
            if (!sscanf(token, "%hu", &swath)) {
              if (line == 1) { // Check if it is the swath header
                if (!(strncmp(token, "\"Swath\"", 7) == 0 || strncmp(token, "\"swath\"", 7) == 0 || strncmp(token, "Swath", 5) == 0 || strncmp(token, "swath", 5) == 0)) {
                  fprintf(stderr, "The third column in the filter list must be swath number; got %s instead", token);
                  return 1;
                }
              }
              else {
                fprintf(stderr, "%s:%d: Not an integer: %s", list_file, line, token);
                return 1;
              }
            }
            break;

          case 4:
            if (!sscanf(token, "%hu", &tile)) {
              if (line == 1) { // Check if it is the tile header
                if (!(strncmp(token, "\"Tile\"", 6) == 0 || strncmp(token, "\"tile\"", 6) == 0 || strncmp(token, "Tile", 4) == 0 || strncmp(token, "tile", 4) == 0)) {
                  fprintf(stderr, "The fourth column in the filter list must be tile number; got %s instead", token);
                  return 1;
                }
              }
              else {
                fprintf(stderr, "%s:%d: Not an integer: %s", list_file, line, token);
                return 1;
              }
            }

          default:
            goto end_of_scan;  // ignore extra columns
        }  // switch (tn)

        token = strtok(NULL, ",");
      }  // while (token)

end_of_scan:
      if (tn < 4) {
        fprintf(stderr, "%s:%d: The filter list must contain four tile co-ordinates: lane, surface, swath, and tile number", list_file, line);
        return 1;
      }

      // Check values
      if (lane == 0 && surface == 0 && swath == 0 && tile == 0) {
        continue;  // Line 1 must be a header line
      }

      if (lane < 1 || lane > 4) {
        fprintf(stderr, "%s:%d: Invalid lane number %hu", list_file, line, lane);
        return 1;
      }
      if (surface < 1 || surface > 2) {
        fprintf(stderr, "%s:%d: Invalid surface number %hu", list_file, line, surface);
        return 1;
      }
      if (swath < 1 || swath > 6) {
        fprintf(stderr, "%s:%d: Invalid swath number %hu", list_file, line, swath);
        return 1;
      }
      if (tile < 1 || tile > 78) {
        fprintf(stderr, "%s:%d: Invalid tile number %hu", list_file, line, tile);
        return 1;
      }

      // All good
      filter[lane - 1][surface - 1][swath - 1][tile - 1] ^= 1;
    }  // read filter list


    if (fread(magic, 1, 2, stdin) == 2) {
      if (*(uint16_t*)magic == 0x8b1f) {
        fprintf(stderr, "Could not read compressed FASTQ from stdin\n");
        return 1;
      }
      if (magic[0] != '@') {
        fprintf(stderr, "Didn't see FASTQ on stdin\n");
        return 1;
      }
      fseek(stdin, 0, SEEK_SET);
    }
    else {
      if (feof(stdin)) {
        fprintf(stderr, "less than 2 bytes written to stdin\n");
      }
      else {
        fprintf(stderr, "Error reading from stdin\n");
      }
      return 1;
    }

    scan_tiles(stdin);
  }
  else {
    for (; optind < argc; optind++) {
      FILE*  fin = fopen(argv[optind], "rb");
      if (fin == NULL) {
        fprintf(stderr, "Could not open %s\n", argv[optind]);
        continue;
      }
      if (fread(magic, 1, 2, fin) == 2) {
        if (*(uint16_t*)magic == 0x8b1f) {
          fprintf(stderr, "Can't read compressed file %s\n", argv[optind]);
          return 1;
        }
        if (magic[0] != '@') {
          fprintf(stderr, "Is %s a FASTQ file?\n", argv[optind]);
          return 1;
        }
        fseek(fin, 0, SEEK_SET);
      }
      else {
        if (feof(fin)) {
          fprintf(stderr, "File `%s` is less than 2 bytes long\n", argv[optind]);
        }
        else {
          fprintf(stderr, "Error reading '%s'\n", argv[optind]);
        }
        return 1;
      }

      scan_tiles(fin);

      fclose(fin);
    }
  }


  return 0;
}



