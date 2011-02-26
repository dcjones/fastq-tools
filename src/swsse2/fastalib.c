/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "swsse2.h"
#include "fastalib.h"

#define READ_BUFFER_SIZE (128 * 1024)
#define SEQ_NAME_SIZE    (128)

FASTA_LIB *openLib (char *file, int pad)
{
    FILE *fp;

    FASTA_LIB *lib;

    if ((fp = fopen (file, "r")) == NULL) {
        fprintf (stderr, "Unable to open file %s\n", file);
        exit (-1);
    }

    lib = (FASTA_LIB *) malloc (sizeof (FASTA_LIB));
    if (!lib) {
        fprintf (stderr, "Unable to allocate memory for library header\n");
        exit (-1);
    }

    lib->readBuffer = (char *) malloc (READ_BUFFER_SIZE);
    if (!lib->readBuffer) {
        fprintf (stderr, "Unable to allocate memory for read buffer\n");
        exit (-1);
    }

    lib->seqBuffer = (unsigned char *) malloc (MAX_SEQ_LENGTH);
    if (!lib->seqBuffer) {
        fprintf (stderr, "Unable to allocate memory for sequence\n");
        exit (-1);
    }

    lib->seqName = (char *) malloc (SEQ_NAME_SIZE);
    if (!lib->seqName) {
        fprintf (stderr, "Unable to allocate memory for sequence name\n");
        exit (-1);
    }

    lib->size = (int) fread (lib->readBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
    if (lib->size == 0 && !feof (fp)) {
        fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
        exit (-1);
    }

    lib->pos = 0;

    lib->fp = fp;

    lib->sequences = 0;
    lib->residues = 0;

    lib->pad;

    return lib;
}

static int
readNextBlock (FASTA_LIB *lib)
{
    FILE *fp = lib->fp;
    size_t size;

    size = fread (lib->readBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
    if (lib->size == 0 && !feof (fp)) {
        fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
        exit (-1);
    }

    lib->pos = 0;
    lib->size = (int) size;

    return lib->size;
}

unsigned char *
nextSeq (FASTA_LIB *lib, int *length)
{
    int inx;
    int size;
    int done;
    int len;

    char *name = lib->seqName;
    unsigned char *seq = lib->seqBuffer;

    /* check if we are at the end of the library */
    if (lib->size == 0) {
        *length = 0;
        return NULL;
    }

    if (lib->pos == lib->size) {
        readNextBlock (lib);
    }

    inx = lib->pos;

    /* check for the start of a sequence */
    if (lib->readBuffer[inx] != '>') {
        fprintf (stderr, "Error parsing fasta file expecting > found %c\n",
            lib->readBuffer[inx]);
        exit (-1);
    }

    ++inx;

    /* read in the sequence name */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {
            size = readNextBlock (lib);
            if (size == 0) {
                *length = 0;
                return NULL;
            }
            inx = lib->pos;
        } else if (lib->readBuffer[inx] == '\n') {
            *name = '\0';
            done = 1;
        } else if (len < SEQ_NAME_SIZE - 1) {
            *name++ = lib->readBuffer[inx];
            len++;
        }
        ++inx;
    } while (!done);

    lib->pos = inx;

    /* read in the sequence */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {
            size = readNextBlock (lib);
            if (size == 0) {
                *seq = '\0';
                done = 1;
            }
            inx = 0;
        } else if (isspace(lib->readBuffer[inx])) {
            ++inx;
        } else if (lib->readBuffer[inx] == '>') {
            *seq = '\0';
            done = 1;
        } else if (len >= MAX_SEQ_LENGTH) {
            fprintf (stderr, "Sequence %s exceeds maximum length\n", 
                lib->seqName);
            exit (-1);
        } else {
            int value = AMINO_ACID_VALUE[lib->readBuffer[inx]];
            if (value == -1) {
                fprintf (stderr, "Unknown amino acid %c in sequence %s\n",
                    lib->readBuffer[inx], lib->seqName);
                exit (-1);
            }
            *seq++ = (char) value;
            inx++;
            len++;
        }
    } while (!done);

    lib->pos = inx;
    *length = len;

    lib->sequences++;
    lib->residues += len;

    /*  check if we need to pad the sequence to a multiple of 16  */
    if (lib->pad) {
        inx = 16 - (len % 16);
        while (inx--) {
            *seq++ = ALPHA_SIZE;
        }
        *seq = '\0';
    }

    return lib->seqBuffer;
}

void closeLib (FASTA_LIB *lib)
{
    fclose (lib->fp);
    
    free (lib->readBuffer);
    free (lib->seqBuffer);
    free (lib->seqName);

    free (lib);
}
