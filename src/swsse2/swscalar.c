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

#include "swsse2.h"
#include "swscalar.h"

typedef struct {
    int  *pMatrix;
    int  *pH;
    int  *pE;
    int  *pData;
} SwScalarData;

int
swScalar (unsigned char   *querySeq,
          int              queryLength,
          unsigned char   *dbSeq,
          int              dbLength,
          unsigned short   gapOpen,
          unsigned short   gapExtend,
          int             *pMatrix,
          int             *pH,
          int             *pE);

void *
swScalarInit(unsigned char   *querySeq,
             int              queryLength,
             signed char     *matrix)
{
    int i, j;
    int nCount;

    int *pi;

    SwScalarData *pSwData;
 
    /* remove unreferenced warnings */
    querySeq;

    pSwData = (SwScalarData *) malloc (sizeof (SwScalarData));
    if (!pSwData) {
        fprintf (stderr, "Unable to allocate memory for SW data\n");
        exit (-1);
    }

    nCount = ALPHA_SIZE * ALPHA_SIZE +        /* matrix */
             (queryLength * 3);               /* vH1, vH2 and vE */

    pSwData->pData = (int *) calloc (nCount, sizeof (int));
    if (!pSwData->pData) {
        fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
        exit (-1);
    }

    pSwData->pMatrix = pSwData->pData;

    pSwData->pH = pSwData->pMatrix + ALPHA_SIZE * ALPHA_SIZE;
    pSwData->pE  = pSwData->pH + queryLength;

    /* Conver the scoring matrix to ints */
    pi = pSwData->pMatrix;
    for (i = 0; i < ALPHA_SIZE; ++i) {
        for (j = 0; j < ALPHA_SIZE; ++j) {
            *pi++ = (int) *matrix++;
        }
    }

    return pSwData;
}

void swScalarScan (unsigned char   *querySeq,
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

    SwScalarData *scalarData = (SwScalarData *) swData;

    dbSeq = nextSeq (dbLib, &dbLen);
    while (dbLen > 0) {

        score = swScalar (querySeq, queryLength, 
                          dbSeq, dbLen, 
                          gapInit, gapExt, 
                          scalarData->pMatrix,
                          scalarData->pH,
                          scalarData->pE);

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
swScalarComplete(void *pSwData)
{
    SwScalarData *pScalarData = (SwScalarData *) pSwData;

    free (pScalarData->pData);
    free (pScalarData);
}

int
swScalar(unsigned char   *querySeq,
         int              queryLength,
         unsigned char   *dbSeq,
         int              dbLength,
         unsigned short   gapOpen,
         unsigned short   gapExtend,
         int             *pMatrix,
         int             *pH,
         int             *pE)
{
    int     i, j;

    int E, F, H;
    int prevH;

    int maxScore = 0;

    int *pScore;

    /* Zero out the storage vector */
    for (i = 0; i < queryLength; i++)
    {
        *(pE + i) = 0;
        *(pH + i) = 0;
    }

    for (i = 0; i < dbLength; ++i)
    {
        /* fetch first data asap. */
        pScore = pMatrix + dbSeq[i] * ALPHA_SIZE;

        /* zero out F. */
        F = 0;

        /* load the next h value */
        prevH = 0;

        for (j = 0; j < queryLength; j++)
        {
            /* load E values */
            E = *(pE + j);

            /* add score to H */
            H = prevH + *(pScore + querySeq[j]);
            if (H < 0) {
                H = 0;
            }

            /* Update highest score encountered this far */
            if (H > maxScore) {
                maxScore = H;
            }

            /* get max from H, E and F */
            if (E > H) {
                H = E;
            }
            if (F > H) {
                H = F;
            }

            /* save H values */
            prevH = *(pH + j);
            *(pH + j) = H;

            H = H - gapOpen;

            /* update E value */
            E = E - gapExtend;
            if (H > E) {
                E = H;
            }
            if (E < 0) {
                E = 0;
            }

            /* update vF value */
            F = F - gapExtend;
            if (H > F) {
                F = H;
            }
            if (F < 0) {
                F = 0;
            }

            /* save vE values */
            *(pE + j) = E;
        }
    }

    /* return largest score */
    return maxScore;
}
