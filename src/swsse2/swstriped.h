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

#ifndef INCLUDE_SWSTRIPED_H
#define INCLUDE_SWSTRIPED_H

#include <emmintrin.h>

#include "swsse2.h"
#include "fastalib.h"

SW_DATA *
swStripedInit (unsigned char   *querySeq,
               int              queryLength,
               signed char     *matrix);


void 
swStripedScan (unsigned char   *querySeq,
               int              queryLength,
               FASTA_LIB       *dbLib,
               void            *swData,
               SEARCH_OPTIONS  *options,
               SCORE_LIST      *scores);

                    
void
swStripedComplete (SW_DATA *pSwData);

int
swStripedWord (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               __m128i         *queryProf,
               __m128i         *pvH1,
               __m128i         *pvH2,
               __m128i         *pvE);


int
swStripedByte (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               __m128i         *queryProf,
               __m128i         *pvH1,
               __m128i         *pvH2,
               __m128i         *pvE,
               unsigned short   bias);

typedef struct {
    __m128i        *pvbQueryProf;
    __m128i        *pvsQueryProf;
    __m128i        *pvH1;
    __m128i        *pvH2;
    __m128i        *pvE;
    unsigned char  *pData;
    unsigned short  bias;
} SwStripedData;



#endif /* INCLUDE_SWSTRIPED_H */
