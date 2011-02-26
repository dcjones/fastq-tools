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

#ifndef INCLUDE_FASTALIB_H
#define INCLUDE_FASTALIB_H

#include <stdio.h>

#define MAX_SEQ_LENGTH (64 * 1024)

typedef struct {
    char *readBuffer;

    char *seqName;
    unsigned char *seqBuffer;

    int pos;
    int size;

    FILE *fp;

    int sequences;
    int residues;

    int pad;
} FASTA_LIB;

FASTA_LIB *openLib (char *file, int pad);
void closeLib (FASTA_LIB *lib);

unsigned char *nextSeq (FASTA_LIB *lib, int *length);

#define seqName(LIB) (LIB->seqName)

#endif /* INCLUDE_FASTALIB_H */
