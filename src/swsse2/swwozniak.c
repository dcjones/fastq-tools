/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/* Written by Michael Farrar, 2006.
   Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

/* Implementation of the Wozniak "vertical" vectorization
   strategy for Smith-Waterman comparison, Wozniak (1997) Comp.
   Appl. Biosci. 13:145-150
*/

#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>

#include "swsse2.h"
#include "swwozniak.h"

#define MATRIX_ROW_SIZE 32
#define MATRIX_SIZE     (MATRIX_ROW_SIZE * (ALPHA_SIZE + 1))

typedef struct {
    char           *pbMatrix;
    short          *psMatrix;
    __m128i        *pvHStore;
    __m128i        *pvEStore;
    unsigned char  *pData;
    unsigned short  bias;
} SwWozniakData;

int
swWozniakWord (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               short           *pMatrix,
               __m128i         *pvHStore,
               __m128i         *pvEStore);

int
swWozniakByte (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               char            *pMatrix,
               __m128i         *pvHStore,
               __m128i         *pvEStore,
               unsigned short   bias);

void *
swWozniakInit(unsigned char   *querySeq,
              int              queryLength,
              signed char     *matrix)
{
    int i, j;

    int nCount;

    int lenQryByte;
    int lenQryShort;

    int bias;

    int weight;

    short *ps;
    char *pc;

    signed char *matrixRow;

    size_t aligned;

    SwWozniakData *pSwData;
 
    lenQryByte = (queryLength + 15) / 16 + 2;
    lenQryShort = (queryLength + 7) / 8 + 2;

    pSwData = (SwWozniakData *) malloc (sizeof (SwWozniakData));
    if (!pSwData) {
        fprintf (stderr, "Unable to allocate memory for SW data\n");
        exit (-1);
    }

    nCount = 64 +                             /* slack bytes */
             4 * (ALPHA_SIZE + 1) +           /* byte matrix */
             4 * (ALPHA_SIZE + 1) +           /* short matrix */
             ((queryLength + 16) * 2);        /* vH and vE */


    pSwData->pData = (unsigned char *) calloc (nCount, sizeof (__m128i));
    if (!pSwData->pData) {
        fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
        exit (-1);
    }

    pSwData->pbMatrix = (char *) pSwData->pData;
    pSwData->psMatrix = (short *) (pSwData->pbMatrix + MATRIX_SIZE);

    /* since we might port this to another platform, lets align the data */
    /* to 16 byte boundries ourselves */
    aligned = (size_t) (pSwData->psMatrix + MATRIX_SIZE);
    aligned = (aligned + 15) & ~(0x0f);

    pSwData->pvHStore = (__m128i *) aligned;
    pSwData->pvEStore = pSwData->pvHStore + queryLength + 16;

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

    pc = pSwData->pbMatrix;
    ps = pSwData->psMatrix;

    for (i = 0; i < ALPHA_SIZE; i++) {
        matrixRow = matrix + i * ALPHA_SIZE;

        for (j = 0; j < ALPHA_SIZE; j++) {
            weight = matrixRow[j];
            *pc++ = weight - bias;
            *ps++ = weight;
        }

        for ( ; j < MATRIX_ROW_SIZE; j++) {
            *pc++ = -bias;
            *ps++ = 0;
        }
    }

    /* add the weights for the NULL rows */
    for (j = 0; j < MATRIX_ROW_SIZE; j++) {
        *pc++ = -bias;
        *ps++ = 0;
    }

    pSwData->bias = (unsigned short) -bias;

    return pSwData;
}

void swWozniakScan (unsigned char   *querySeq,
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

    SwWozniakData *wozniakData = (SwWozniakData *) swData;

    dbSeq = nextSeq (dbLib, &dbLen);
    while (dbLen > 0) {

        score = swWozniakByte (querySeq, queryLength, 
                               dbSeq, dbLen, 
                               gapInit, gapExt, 
                               wozniakData->pbMatrix,
                               wozniakData->pvHStore,
                               wozniakData->pvEStore,
                               wozniakData->bias);

        /* check if needs a run with higher precision */
        if (score >= 255) {
            score = swWozniakWord (querySeq, queryLength, 
                                   dbSeq, dbLen, 
                                   gapInit, gapExt, 
                                   wozniakData->psMatrix,
                                   wozniakData->pvHStore,
                                   wozniakData->pvEStore);
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
swWozniakComplete(void *pSwData)
{
    SwWozniakData *pWozniakData = (SwWozniakData *) pSwData;

    free (pWozniakData->pData);
    free (pWozniakData);
}


int
swWozniakWord (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               short           *pMatrix,
               __m128i         *pvHStore,
               __m128i         *pvEStore)
{
    int      i, j, k, l, m;
    int      score;

    short   *pScore;

    __m128i  vE, vF, vH;
    __m128i  vEUp, vHUp1, vHUp2;

    __m128i  vMaxScore;
    __m128i  vGapOpen;
    __m128i  vGapExtend;
    __m128i  vScore;

    __m128i  vMin;
    __m128i  vMinimums;
    __m128i  vTemp;

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

    for (i = 0; i < queryLength + 8; i++)
    {
        _mm_store_si128 (pvEStore + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }

    pScore = (short *) &vScore;

    for (i = 0; i < dbLength; i += 8, dbSeq += 8)
    {
        /* zero lots of stuff. */
        vE     = _mm_shuffle_epi32 (vMinimums, 0);
        vF     = _mm_shuffle_epi32 (vMinimums, 0);
        vH     = _mm_shuffle_epi32 (vMinimums, 0);
        vHUp2  = _mm_shuffle_epi32 (vMinimums, 0);

        vScore = _mm_xor_si128 (vScore, vScore);

        for (j = 0; j < 8; ++j)
        {
            for (k = 0; k <= j; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }
            for ( ; k < 8; ++k) {
                pScore[k] = 0;
            }

            /* load values of vE and vH from previous row (one unit up) */
            vEUp    = _mm_load_si128 (pvEStore + j);
            vHUp1   = _mm_load_si128 (pvHStore + j);

            /* shift into place so we have complete vE and vH vectors */
            /* that refer to the values one unit up from each cell */
            /* that we are currently working on. */
            vTemp   = _mm_slli_si128 (vE, 2);
            vEUp    = _mm_srli_si128 (vEUp, 14);
            vEUp    = _mm_or_si128 (vEUp, vTemp);

            vTemp   = _mm_slli_si128 (vH, 2);
            vHUp1   = _mm_srli_si128 (vHUp1, 14);
            vHUp1   = _mm_or_si128 (vHUp1, vTemp);

            /* do the dynamic programming */

            /* update vE value */
            vE    = _mm_subs_epi16 (vEUp, vGapExtend);
            vTemp = _mm_subs_epi16 (vHUp1, vGapOpen);
            vE    = _mm_max_epi16 (vE, vTemp);

            /* update vF value */
            vF    = _mm_subs_epi16 (vF, vGapExtend);
            vTemp = _mm_subs_epi16 (vH, vGapOpen);
            vF    = _mm_max_epi16  (vF, vTemp);

            /* add score to vH */
            vH = _mm_adds_epi16 (vHUp2, vScore);

            /* set vH to max of vH, vE, vF */
            vH = _mm_max_epi16 (vH, vE);
            vH = _mm_max_epi16 (vH, vF);

            /* Save value to use for next diagonal vH */
            vHUp2 = vHUp1;

            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epi16 (vMaxScore, vH);
        }

        for (l = 0; j < queryLength; ++j, ++l)
        {
            for (k = 0; k < 8; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }

            /* load values of vE and vH from previous row (one unit up) */
            vEUp    = _mm_load_si128 (pvEStore + j);
            vHUp1   = _mm_load_si128 (pvHStore + j);

            /* save old values of vE and vH to use on next row */
            _mm_store_si128 (pvEStore + l, vE);
            _mm_store_si128 (pvHStore + l, vH);
            
            /* shift into place so we have complete vE and vH vectors */
            /* that refer to the values one unit up from each cell */
            /* that we are currently working on. */
            vTemp = _mm_slli_si128 (vE, 2);
            vEUp  = _mm_srli_si128 (vEUp, 14);
            vEUp  = _mm_or_si128 (vEUp, vTemp);

            vTemp = _mm_slli_si128 (vH, 2);
            vHUp1 = _mm_srli_si128 (vHUp1, 14);
            vHUp1 = _mm_or_si128 (vHUp1, vTemp);

            /* do the dynamic programming */

            /* update vE value */
            vE    = _mm_subs_epi16 (vEUp, vGapExtend);
            vTemp = _mm_subs_epi16 (vHUp1, vGapOpen);
            vE    = _mm_max_epi16 (vE, vTemp);

            /* update vF value */
            vF    = _mm_subs_epi16 (vF, vGapExtend);
            vTemp = _mm_subs_epi16 (vH, vGapOpen);
            vF    = _mm_max_epi16 (vF, vTemp);

            /* add score to vH */
            vH = _mm_adds_epi16(vHUp2, vScore);

            /* set vH to max of vH, vE, vF */
            vH = _mm_max_epi16 (vH, vE);
            vH = _mm_max_epi16 (vH, vF);

            /* Save value to use for next diagonal vH */
            vHUp2 = vHUp1;

            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epi16 (vMaxScore, vH);
        }

        for (m = 0 ; m < 7; ++j, ++l, ++m)
        {
            for (k = 0; k <= m; ++k) {
                pScore[k] = 0;
            }
            for ( ; k < 8; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }

            /* save old values of vE and vH to use on next row */
            _mm_store_si128 (pvEStore + l, vE);
            _mm_store_si128 (pvHStore + l, vH);

            /* v_score_load contains all zeros */
            vTemp = _mm_slli_si128 (vE, 2);
            vEUp  = _mm_or_si128 (vMin, vTemp);
            vTemp = _mm_slli_si128 (vH, 2);
            vHUp1 = _mm_or_si128 (vMin, vTemp);

            /* do the dynamic programming */

            /* update vE value */
            vE    = _mm_subs_epi16 (vEUp, vGapExtend);
            vTemp = _mm_subs_epi16 (vHUp1, vGapOpen);
            vE    = _mm_max_epi16 (vE, vTemp);

            /* update vF value */
            vF    = _mm_subs_epi16 (vF, vGapExtend);
            vTemp = _mm_subs_epi16 (vH, vGapOpen);
            vF    = _mm_max_epi16 (vF, vTemp);

            /* add score to vH */
            vH = _mm_adds_epi16 (vHUp2, vScore);

            /* set vH to max of vH, vE, vF */
            vH = _mm_max_epi16 (vH, vE);
            vH = _mm_max_epi16 (vH, vF);

            /* Save value to use for next diagonal vH */
            vHUp2 = vHUp1;

            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epi16(vMaxScore,vH);
        }

        _mm_store_si128 (pvEStore + l, vE);
        _mm_store_si128 (pvHStore + l, vH);
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
swWozniakByte (unsigned char   *querySeq,
               int              queryLength,
               unsigned char   *dbSeq,
               int              dbLength,
               unsigned short   gapOpen,
               unsigned short   gapExtend,
               char            *pMatrix,
               __m128i         *pvHStore,
               __m128i         *pvEStore,
               unsigned short   bias)
{
    int      i, j, k, l, m;
    int      score;

    int      dup;

    char    *pScore;

    __m128i  vE, vF, vH;
    __m128i  vEUp, vHUp1, vHUp2;

    __m128i  vMaxScore;
    __m128i  vBias;
    __m128i  vGapOpen;
    __m128i  vGapExtend;
    __m128i  vScore;

    __m128i  vTemp;

    /* remove unreferenced warning */
    querySeq;

    /* Load the bias to all elements of a constant */
    dup   = (bias << 8) | (bias & 0x00ff);
    vBias = _mm_insert_epi16 (vBias, dup, 0);
    vBias = _mm_shufflelo_epi16 (vBias, 0);
    vBias = _mm_shuffle_epi32 (vBias, 0);

    /* Load gap opening penalty to all elements of a constant */
    dup      = (gapOpen << 8) | (gapOpen & 0x00ff);
    vGapOpen = _mm_insert_epi16 (vGapOpen, dup, 0);
    vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
    vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

    /* Load gap extension penalty to all elements of a constant */
    dup        = (gapExtend << 8) | (gapExtend & 0x00ff);
    vGapExtend = _mm_insert_epi16 (vGapExtend, dup, 0);
    vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
    vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

    vScore = _mm_xor_si128 (vScore, vScore);
    vMaxScore = _mm_xor_si128 (vMaxScore, vMaxScore);

    for (i = 0; i < queryLength + 16; i++)
    {
        _mm_store_si128 (pvEStore + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }

    pScore = (char *) &vScore;

    for (i = 0; i < dbLength; i += 16, dbSeq += 16)
    {
        // zero lots of stuff.
        vE     = _mm_xor_si128 (vE, vE);
        vF     = _mm_xor_si128 (vF, vF);
        vH     = _mm_xor_si128 (vH, vH);
        vHUp2  = _mm_xor_si128 (vHUp2, vHUp2);

        vScore = _mm_xor_si128 (vScore, vScore);

        for (j = 0; j < 16; ++j)
        {
            for (k = 0; k <= j; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }
            for ( ; k < 16; ++k) {
                pScore[k] = (char) bias;
            }

            // load values of vE and vH from previous row (one unit up)
            vEUp    = _mm_load_si128 (pvEStore + j);
            vHUp1   = _mm_load_si128 (pvHStore + j);

            // shift into place so we have complete vE and vH vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            vTemp   = _mm_slli_si128 (vE, 1);
            vEUp    = _mm_srli_si128 (vEUp, 15);
            vEUp    = _mm_or_si128 (vEUp, vTemp);

            vTemp   = _mm_slli_si128 (vH, 1);
            vHUp1   = _mm_srli_si128 (vHUp1, 15);
            vHUp1   = _mm_or_si128 (vHUp1, vTemp);

            // do the dynamic programming

            // update vE value
            vE    = _mm_subs_epu8 (vEUp, vGapExtend);
            vTemp = _mm_subs_epu8 (vHUp1, vGapOpen);
            vE    = _mm_max_epu8 (vE, vTemp);

            // update vF value
            vF    = _mm_subs_epu8 (vF, vGapExtend);
            vTemp = _mm_subs_epu8 (vH, vGapOpen);
            vF    = _mm_max_epu8  (vF, vTemp);

            // add score to vH
            vH = _mm_adds_epu8 (vHUp2, *((__m128i *) pScore));
            vH = _mm_subs_epu8 (vH, vBias);

            // set vH to max of vH, vE, vF
            vH = _mm_max_epu8 (vH, vE);
            vH = _mm_max_epu8 (vH, vF);

            // Save value to use for next diagonal vH
            vHUp2 = vHUp1;

            // Update highest score encountered this far
            vMaxScore = _mm_max_epu8 (vMaxScore, vH);
        }

        for (l = 0; j < queryLength; ++j, ++l)
        {
            for (k = 0; k < 16; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }

            // load values of vE and vH from previous row (one unit up)
            vEUp    = _mm_load_si128 (pvEStore + j);
            vHUp1   = _mm_load_si128 (pvHStore + j);

            // save old values of vE and vH to use on next row
            _mm_store_si128 (pvEStore + l, vE);
            _mm_store_si128 (pvHStore + l, vH);
            
            // shift into place so we have complete vE and vH vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            vTemp = _mm_slli_si128 (vE, 1);
            vEUp  = _mm_srli_si128 (vEUp, 15);
            vEUp  = _mm_or_si128 (vEUp, vTemp);

            vTemp = _mm_slli_si128 (vH, 1);
            vHUp1 = _mm_srli_si128 (vHUp1, 15);
            vHUp1 = _mm_or_si128 (vHUp1, vTemp);

            // do the dynamic programming

            // update vE value
            vE    = _mm_subs_epu8 (vEUp, vGapExtend);
            vTemp = _mm_subs_epu8 (vHUp1, vGapOpen);
            vE    = _mm_max_epu8 (vE, vTemp);

            // update vF value
            vF    = _mm_subs_epu8 (vF, vGapExtend);
            vTemp = _mm_subs_epu8 (vH, vGapOpen);
            vF    = _mm_max_epu8 (vF, vTemp);

            // add score to vH
            vH = _mm_adds_epu8(vHUp2, vScore);
            vH = _mm_subs_epu8 (vH, vBias);

            // set vH to max of vH, vE, vF
            vH = _mm_max_epu8 (vH, vE);
            vH = _mm_max_epu8 (vH, vF);

            // Save value to use for next diagonal vH
            vHUp2 = vHUp1;

            // Update highest score encountered this far
            vMaxScore = _mm_max_epu8 (vMaxScore, vH);
        }

        for (m = 0 ; m < 15; ++j, ++l, ++m)
        {
            for (k = 0; k <= m; ++k) {
                pScore[k] = (char) bias;
            }
            for ( ; k < 16; ++k) {
                int matrixOffset = *(dbSeq + k) * MATRIX_ROW_SIZE;
                pScore[k] = *(pMatrix + matrixOffset + *(querySeq + j - k));
            }

            // save old values of vE and vH to use on next row
            _mm_store_si128 (pvEStore + l, vE);
            _mm_store_si128 (pvHStore + l, vH);

            // v_score_load contains all zeros
            vEUp  = _mm_slli_si128 (vE, 1);
            vHUp1 = _mm_slli_si128 (vH, 1);

            // do the dynamic programming

            // update vE value
            vE    = _mm_subs_epu8 (vEUp, vGapExtend);
            vTemp = _mm_subs_epu8 (vHUp1, vGapOpen);
            vE    = _mm_max_epu8 (vE, vTemp);

            // update vF value
            vF    = _mm_subs_epu8 (vF, vGapExtend);
            vTemp = _mm_subs_epu8 (vH, vGapOpen);
            vF    = _mm_max_epu8 (vF, vTemp);

            // add score to vH
            vH = _mm_adds_epu8 (vHUp2, vScore);
            vH = _mm_subs_epu8 (vH, vBias);

            // set vH to max of vH, vE, vF
            vH = _mm_max_epu8 (vH, vE);
            vH = _mm_max_epu8 (vH, vF);

            // Save value to use for next diagonal vH
            vHUp2 = vHUp1;

            // Update highest score encountered this far
            vMaxScore = _mm_max_epu8(vMaxScore,vH);
        }

        _mm_store_si128 (pvEStore + l, vE);
        _mm_store_si128 (pvHStore + l, vH);
    }

    // find largest score in the vMaxScore vector
    vTemp = _mm_srli_si128 (vMaxScore, 8);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 4);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 2);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 1);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);

    // store in temporary variable
    score = (short) _mm_extract_epi16 (vMaxScore, 0);
    score = score & 0x00ff;

    //  check if we might have overflowed
    if (score + bias >= 255)
    {
        score = 255;
    }


    // return largest score
    return score;
}
