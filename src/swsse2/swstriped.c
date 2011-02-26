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
#include <emmintrin.h>

#include "swsse2.h"
#include "swstriped.h"

void *
swStripedInit(unsigned char   *querySeq,
              int              queryLength,
              signed char     *matrix)
{
    int i, j, k;

    int segSize;
    int nCount;

    int bias;

    int lenQryByte;
    int lenQryShort;

    int weight;

    short *ps;
    char *pc;

    signed char *matrixRow;

    size_t aligned;

    SwStripedData *pSwData;
 
    lenQryByte = (queryLength + 15) / 16;
    lenQryShort = (queryLength + 7) / 8;

    pSwData = (SwStripedData *) malloc (sizeof (SwStripedData));
    if (!pSwData) {
        fprintf (stderr, "Unable to allocate memory for SW data\n");
        exit (-1);
    }

    nCount = 64 +                             /* slack bytes */
             lenQryByte * ALPHA_SIZE +        /* query profile byte */
             lenQryShort * ALPHA_SIZE +       /* query profile short */
             (lenQryShort * 3);               /* vH1, vH2 and vE */

    pSwData->pData = (unsigned char *) calloc (nCount, sizeof (__m128i));
    if (!pSwData->pData) {
        fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
        exit (-1);
    }

    /* since we might port this to another platform, lets align the data */
    /* to 16 byte boundries ourselves */
    aligned = ((size_t) pSwData->pData + 15) & ~(0x0f);

    pSwData->pvbQueryProf = (__m128i *) aligned;
    pSwData->pvsQueryProf = pSwData->pvbQueryProf + lenQryByte * ALPHA_SIZE;

    pSwData->pvH1 = pSwData->pvsQueryProf + lenQryShort * ALPHA_SIZE;
    pSwData->pvH2 = pSwData->pvH1 + lenQryShort;
    pSwData->pvE  = pSwData->pvH2 + lenQryShort;

    /* Use a scoring profile for the SSE2 implementation, but the layout
     * is a bit strange.  The scoring profile is parallel to the query, but is
     * accessed in a stripped pattern.  The query is divided into equal length
     * segments.  The number of segments is equal to the number of elements
     * processed in the SSE2 register.  For 8-bit calculations, the query will
     * be divided into 16 equal length parts.  If the query is not long enough
     * to fill the last segment, it will be filled with neutral weights.  The
     * first element in the SSE register will hold a value from the first segment,
     * the second element of the SSE register will hold a value from the
     * second segment and so on.  So if the query length is 288, then each
     * segment will have a length of 18.  So the first 16 bytes will  have
     * the following weights: Q1, Q19, Q37, ... Q271; the next 16 bytes will
     * have the following weights: Q2, Q20, Q38, ... Q272; and so on until
     * all parts of all segments have been written.  The last seqment will
     * have the following weights: Q18, Q36, Q54, ... Q288.  This will be
     * done for the entire alphabet.
     */

    /* Find the bias to use in the substitution matrix */
    bias = 127;
    for (i = 0; i < ALPHA_SIZE * ALPHA_SIZE; i++) {
        if (matrix[i] < bias) {
            bias = matrix[i];
        }
    }
    if (bias > 0) {
        bias = 0;
    }

    /* Fill in the byte query profile */
    pc = (char *) pSwData->pvbQueryProf;
    segSize = (queryLength + 15) / 16;
    nCount = segSize * 16;
    for (i = 0; i < ALPHA_SIZE; ++i) {
        matrixRow = matrix + i * ALPHA_SIZE;
        for (j = 0; j < segSize; ++j) {
            for (k = j; k < nCount; k += segSize) {
                if (k >= queryLength) {
                    weight = 0;
                } else {
                    weight = matrixRow[*(querySeq + k)];
                }
                *pc++ = (char) (weight - bias);
            }
        }
    }

    /* Fill in the short query profile */
    ps = (short *) pSwData->pvsQueryProf;
    segSize = (queryLength + 7) / 8;
    nCount = segSize * 8;
    for (i = 0; i < ALPHA_SIZE; ++i) {
        matrixRow = matrix + i * ALPHA_SIZE;
        for (j = 0; j < segSize; ++j) {
            for (k = j; k < nCount; k += segSize) {
                if (k >= queryLength) {
                    weight = 0;
                } else {
                    weight = matrixRow[*(querySeq + k)];
                }
                *ps++ = (unsigned short) weight;
            }
        }
    }

    pSwData->bias = (unsigned short) -bias;

    return pSwData;
}

void swStripedScan (unsigned char   *querySeq,
                    int              queryLength,
                    FASTA_LIB       *dbLib,
                    void            *swData,
                    SEARCH_OPTIONS  *options,
                    SCORE_LIST      *scores)
{
    int score;

    int threshold = options->threshold;

    unsigned char *dbSeq;
    int dbLen;

    int gapInit = -(options->gapInit + options->gapExt);
    int gapExt  = -options->gapExt;

    SwStripedData *stripedData = (SwStripedData *) swData;

    dbSeq = nextSeq (dbLib, &dbLen);
    while (dbLen > 0) {

        score = swStripedByte (querySeq, queryLength, 
                               dbSeq, dbLen, 
                               gapInit, gapExt, 
                               stripedData->pvbQueryProf,
                               stripedData->pvH1,
                               stripedData->pvH2,
                               stripedData->pvE,
                               stripedData->bias);

        /* check if needs a run with higher precision */
        if (score >= 255) {
            score = swStripedWord (querySeq, queryLength, 
                                   dbSeq, dbLen, 
                                   gapInit, gapExt, 
                                   stripedData->pvsQueryProf,
                                   stripedData->pvH1,
                                   stripedData->pvH2,
                                   stripedData->pvE);
        }

        if (score >= threshold) {
            int minScore = insertList (scores, score, seqName (dbLib));
            if (minScore >= threshold) {
                threshold = minScore;
            }
        }

        dbSeq = nextSeq (dbLib, &dbLen);
    }
}

void
swStripedComplete(void *pSwData)
{
    SwStripedData *pStripedData = (SwStripedData *) pSwData;

    free (pStripedData->pData);
    free (pStripedData);
}

int
swStripedWord(unsigned char   *querySeq,
              int              queryLength,
              unsigned char   *dbSeq,
              int              dbLength,
              unsigned short   gapOpen,
              unsigned short   gapExtend,
              __m128i         *pvQueryProf,
              __m128i         *pvHLoad,
              __m128i         *pvHStore,
              __m128i         *pvE)
{
    int     i, j;
    int     score;

    int     cmp;
    int     iter = (queryLength + 7) / 8;
    
    __m128i *pv;

    __m128i vE, vF, vH;

    __m128i vMaxScore;
    __m128i vGapOpen;
    __m128i vGapExtend;

    __m128i vMin;
    __m128i vMinimums;
    __m128i vTemp;

    __m128i *pvScore;

    /* remove unreferenced warning */
    querySeq;

    /* Load gap opening penalty to all elements of a constant */
    vGapOpen = _mm_insert_epi16 (vGapOpen, gapOpen, 0);
    vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
    vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

    /* Load gap extension penalty to all elements of a constant */
    vGapExtend = _mm_insert_epi16 (vGapExtend, gapExtend, 0);
    vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
    vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

    /*  load vMaxScore with the zeros.  since we are using signed */
    /*  math, we will bias the maxscore to -32768 so we have the */
    /*  full range of the short. */
    vMaxScore = _mm_cmpeq_epi16 (vMaxScore, vMaxScore);
    vMaxScore = _mm_slli_epi16 (vMaxScore, 15);

    vMinimums = _mm_shuffle_epi32 (vMaxScore, 0);

    vMin = _mm_shuffle_epi32 (vMaxScore, 0);
    vMin = _mm_srli_si128 (vMin, 14);

    /* Zero out the storage vector */
    for (i = 0; i < iter; i++)
    {
        _mm_store_si128 (pvE + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }

    for (i = 0; i < dbLength; ++i)
    {
        /* fetch first data asap. */
        pvScore = pvQueryProf + dbSeq[i] * iter;

        /* zero out F. */
        vF = _mm_cmpeq_epi16 (vF, vF);
        vF = _mm_slli_epi16 (vF, 15);

        /* load the next h value */
        vH = _mm_load_si128 (pvHStore + iter - 1);
        vH = _mm_slli_si128 (vH, 2);
        vH = _mm_or_si128 (vH, vMin);

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        for (j = 0; j < iter; j++)
        {
            /* load values of vF and vH from previous row (one unit up) */
            vE = _mm_load_si128 (pvE + j);

            /* add score to vH */
            vH = _mm_adds_epi16 (vH, *pvScore++);

            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epi16 (vMaxScore, vH);

            /* get max from vH, vE and vF */
            vH = _mm_max_epi16 (vH, vE);
            vH = _mm_max_epi16 (vH, vF);

            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);

            /* update vE value */
            vH = _mm_subs_epi16 (vH, vGapOpen);
            vE = _mm_subs_epi16 (vE, vGapExtend);
            vE = _mm_max_epi16 (vE, vH);

            /* update vF value */
            vF = _mm_subs_epi16 (vF, vGapExtend);
            vF = _mm_max_epi16 (vF, vH);

            /* save vE values */
            _mm_store_si128 (pvE + j, vE);

            /* load the next h value */
            vH = _mm_load_si128 (pvHLoad + j);
        }

        /* reset pointers to the start of the saved data */
        j = 0;
        vH = _mm_load_si128 (pvHStore + j);

        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = _mm_slli_si128 (vF, 2);
        vF = _mm_or_si128 (vF, vMin);
        vTemp = _mm_subs_epi16 (vH, vGapOpen);
        vTemp = _mm_cmpgt_epi16 (vF, vTemp);
        cmp  = _mm_movemask_epi8 (vTemp);
        while (cmp != 0x0000) 
        {
            vE = _mm_load_si128 (pvE + j);

            vH = _mm_max_epi16 (vH, vF);

            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);

            /*  update vE incase the new vH value would change it */
            vH = _mm_subs_epi16 (vH, vGapOpen);
            vE = _mm_max_epi16 (vE, vH);
            _mm_store_si128 (pvE + j, vE);

            /* update vF value */
            vF = _mm_subs_epi16 (vF, vGapExtend);

            j++;
            if (j >= iter)
            {
                j = 0;
                vF = _mm_slli_si128 (vF, 2);
                vF = _mm_or_si128 (vF, vMin);
            }

            vH = _mm_load_si128 (pvHStore + j);

            vTemp = _mm_subs_epi16 (vH, vGapOpen);
            vTemp = _mm_cmpgt_epi16 (vF, vTemp);
            cmp  = _mm_movemask_epi8 (vTemp);
        }
    }

    /* find largest score in the vMaxScore vector */
    vTemp = _mm_srli_si128 (vMaxScore, 8);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 4);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 2);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);

    /* store in temporary variable */
    score = (short) _mm_extract_epi16 (vMaxScore, 0);

    /* return largest score */
    return score + SHORT_BIAS;
}




int
swStripedByte(unsigned char   *querySeq,
              int              queryLength,
              unsigned char   *dbSeq,
              int              dbLength,
              unsigned short   gapOpen,
              unsigned short   gapExtend,
              __m128i         *pvQueryProf,
              __m128i         *pvHLoad,
              __m128i         *pvHStore,
              __m128i         *pvE,
              unsigned short   bias)
{
    int     i, j;
    int     score;

    int     dup;
    int     cmp;
    int     iter = (queryLength + 15) / 16;
    
    __m128i *pv;

    __m128i vE, vF, vH;

    __m128i vMaxScore;
    __m128i vBias;
    __m128i vGapOpen;
    __m128i vGapExtend;

    __m128i vTemp;
    __m128i vZero;

    __m128i *pvScore;

    /* remove unreferenced warning */
    querySeq;

    /* Load the bias to all elements of a constant */
    dup    = (bias << 8) | (bias & 0x00ff);
    vBias = _mm_insert_epi16 (vBias, dup, 0);
    vBias = _mm_shufflelo_epi16 (vBias, 0);
    vBias = _mm_shuffle_epi32 (vBias, 0);

    /* Load gap opening penalty to all elements of a constant */
    dup    = (gapOpen << 8) | (gapOpen & 0x00ff);
    vGapOpen = _mm_insert_epi16 (vGapOpen, dup, 0);
    vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
    vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

    /* Load gap extension penalty to all elements of a constant */
    dup    = (gapExtend << 8) | (gapExtend & 0x00ff);
    vGapExtend = _mm_insert_epi16 (vGapExtend, dup, 0);
    vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
    vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

    vMaxScore = _mm_xor_si128 (vMaxScore, vMaxScore);

    vZero = _mm_xor_si128 (vZero, vZero);

    /* Zero out the storage vector */
    for (i = 0; i < iter; i++)
    {
        _mm_store_si128 (pvE + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }

    for (i = 0; i < dbLength; ++i)
    {
        /* fetch first data asap. */
        pvScore = pvQueryProf + dbSeq[i] * iter;

        /* zero out F. */
        vF = _mm_xor_si128 (vF, vF);

        /* load the next h value */
        vH = _mm_load_si128 (pvHStore + iter - 1);
        vH = _mm_slli_si128 (vH, 1);

        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        for (j = 0; j < iter; j++)
        {
            /* load values of vF and vH from previous row (one unit up) */
            vE = _mm_load_si128 (pvE + j);

            /* add score to vH */
            vH = _mm_adds_epu8 (vH, *(pvScore + j));
            vH = _mm_subs_epu8 (vH, vBias);

            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epu8 (vMaxScore, vH);

            /* get max from vH, vE and vF */
            vH = _mm_max_epu8 (vH, vE);
            vH = _mm_max_epu8 (vH, vF);

            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);

            /* update vE value */
            vH = _mm_subs_epu8 (vH, vGapOpen);
            vE = _mm_subs_epu8 (vE, vGapExtend);
            vE = _mm_max_epu8 (vE, vH);

            /* update vF value */
            vF = _mm_subs_epu8 (vF, vGapExtend);
            vF = _mm_max_epu8 (vF, vH);

            /* save vE values */
            _mm_store_si128 (pvE + j, vE);

            /* load the next h value */
            vH = _mm_load_si128 (pvHLoad + j);
        }

        /* reset pointers to the start of the saved data */
        j = 0;
        vH = _mm_load_si128 (pvHStore + j);

        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = _mm_slli_si128 (vF, 1);
        vTemp = _mm_subs_epu8 (vH, vGapOpen);
        vTemp = _mm_subs_epu8 (vF, vTemp);
        vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
        cmp  = _mm_movemask_epi8 (vTemp);

        while (cmp != 0xffff) 
        {
            vE = _mm_load_si128 (pvE + j);

            vH = _mm_max_epu8 (vH, vF);

            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);

            /*  update vE incase the new vH value would change it */
            vH = _mm_subs_epu8 (vH, vGapOpen);
            vE = _mm_max_epu8 (vE, vH);
            _mm_store_si128 (pvE + j, vE);

            /* update vF value */
            vF = _mm_subs_epu8 (vF, vGapExtend);

            j++;
            if (j >= iter)
            {
                j = 0;
                vF = _mm_slli_si128 (vF, 1);
            }

            vH = _mm_load_si128 (pvHStore + j);

            vTemp = _mm_subs_epu8 (vH, vGapOpen);
            vTemp = _mm_subs_epu8 (vF, vTemp);
            vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
            cmp  = _mm_movemask_epi8 (vTemp);
        }
    }

    /* find largest score in the vMaxScore vector */
    vTemp = _mm_srli_si128 (vMaxScore, 8);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 4);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 2);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 1);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);

    /* store in temporary variable */
    score = _mm_extract_epi16 (vMaxScore, 0);
    score = score & 0x00ff;

    /*  check if we might have overflowed */
    if (score + bias >= 255)
    {
        score = 255;
    }

    /* return largest score */
    return score;
}
