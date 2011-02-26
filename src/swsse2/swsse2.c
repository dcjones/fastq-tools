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
#include <string.h>

#include <time.h>
#include <sys/timeb.h>

#include "swsse2.h"
#include "matrix.h"
#include "fastalib.h"

#include "swscalar.h"
#include "swwozniak.h"
#ifdef WITH_ROGNES
#include "swrognes.h"
#endif
#include "swstriped.h"

typedef enum { SCALAR, 
               WOZNIAK, 
#ifdef WITH_ROGNES
               ROGNES, 
#endif
               STRIPED 
} SW_TYPES;

const char *SW_IMPLEMENATION[] = {
    "Non-optimized",
    "Wozniak",
#ifdef WITH_ROGNES
    "Rognes",
#endif
    "Striped",
};

typedef struct {
    SW_DATA *(*init) (unsigned char   *querySeq,
                      int              queryLength,
                      signed char     *matrix);
    void     (*scan) (unsigned char   *querySeq,
                      int              queryLength,
                      FASTA_LIB       *dbLib,
                      void            *swData,
                      SEARCH_OPTIONS  *options,
                      SCORE_LIST      *scores);
    void     (*done) (SW_DATA         *pSwData);
} SW_FUNCT_DEFS;

SW_FUNCT_DEFS swFuncs[] = {
    {
        swScalarInit,
        swScalarScan,
        swScalarComplete,
    },
    {
        swWozniakInit,
        swWozniakScan,
        swWozniakComplete,
    },
#ifdef WITH_ROGNES
    {
        swRognesInit,
        swRognesScan,
        swRognesComplete,
    },
#endif
    {
        swStripedInit,
        swStripedScan,
        swStripedComplete,
    },
};

const char AMINO_ACIDS[ALPHA_SIZE] = {
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 
    'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 
    'S', 'T', 'V', 'W', 'X', 'Y', 'Z'
};

const int AMINO_ACID_VALUE[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0,  1,  2,  3,  4,  5,  6,  7,  8, -1,  9, 10, 11, 12, -1,
    13, 14, 15, 16, 17, -1, 18, 19, 20, 21, 22, -1, -1, -1, -1, -1,
    -1,  0,  1,  2,  3,  4,  5,  6,  7,  8, -1,  9, 10, 11, 12, -1,
    13, 14, 15, 16, 17, -1, 18, 19, 20, 21, 22, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

SCORE_LIST *initList (int count);
void freeList (SCORE_LIST *list);

void printResults (SCORE_LIST *list);

void printUsage (void)
{
    printf ("Usage: swsse2 [-h] [-(n|w|r|s)] [-i num] [-e num] [-t num] "
        "[-c num] matrix query db\n");
    printf ("    -h       : this help message\n");
    printf ("    -n       : run a non-vectorized Smith-Waterman search\n");
    printf ("    -w       : run a vectorized Wozniak search\n");
#ifdef WITH_ROGNES
    printf ("    -r       : run a vectorized Rognes search\n");
#else
    printf ("    -r       : run a vectorized Rognes search (NOT SUPPORTED)\n");
#endif
    printf ("    -s       : run a vectorized striped search (default)\n");
    printf ("    -i num   : gap init penalty (default -10)\n");
    printf ("    -e num   : gap extension penalty (default -2)\n");
    printf ("    -t num   : minimum score threshold (default 20)\n");
    printf ("    -c num   : number of scores to be displayed (default 250)\n");
    printf ("    matrix   : scoring matrix file\n");
    printf ("    query    : query sequence file (fasta format)\n");
    printf ("    db       : sequence database file (fasta format)\n");
}

#if 0
int main (int argc, char **argv)
{
    int i;

    int penalty;

    int rptCount = 250;

    SW_TYPES swType = STRIPED;

    char *dbFile = NULL;
    char *queryFile = NULL;
    char *matrixFile = NULL;

    signed char *matrix;

    unsigned char *querySeq;
    int queryLen;

    SCORE_LIST *list;

    FASTA_LIB *queryLib;
    FASTA_LIB *dbLib;

    void *swData;

    struct _timeb startTime;
    struct _timeb endTime;

    double startTimeSec;
    double endTimeSec;

    SEARCH_OPTIONS options;

    if (argc < 4) {
        printUsage ();
        exit (-1);
    }

    /* set the default options */
    options.gapInit = -10;
    options.gapExt = -2;
    options.threshold = 20;

    i = 1;
    while (i < argc) {
        if (i + 3 == argc) {
            /* should be matrix file name */
            matrixFile = argv[i];

        } else if (i + 2 == argc) {
            /* should be query file name */
            queryFile = argv[i];

        } else if (i + 1 == argc) {
            /* should be matrix file name */
            dbFile = argv[i];

        } else {
            /* process arguements */
            switch (argv[i][1]) {
               case 'h':
                   printUsage ();
                   break;
               case 'n':
                   swType = SCALAR;
                   break;
               case 'r':
#ifdef WITH_ROGNES
                   swType = ROGNES;
#else
                   fprintf (stderr, "The ROGNES implementation is not supported\n");
                   exit (-1);
#endif
                   break;
               case 'w':
                   swType = WOZNIAK;
                   break;
               case 's':
                   swType = STRIPED;
                   break;
               case 'i':
                   penalty = atoi (argv[++i]);
                   if (penalty > 0 || penalty < -128) {
                        fprintf (stderr, "Invalid gap init %d\n", penalty);
                        fprintf (stderr, "Valid range is 0 - -128\n");
                        exit (-1);
                   }
                   options.gapInit = (unsigned short) penalty;
                   break;
               case 'e':
                   penalty = atoi (argv[++i]);
                   if (penalty > 0 || penalty < -128) {
                        fprintf (stderr, "Invalid gap extension %d\n", penalty);
                        fprintf (stderr, "Valid range is 0 - -128\n");
                        exit (-1);
                   }
                   options.gapExt = (unsigned short) penalty;
                   break;
               case 't':
                   options.threshold = atoi (argv[++i]);
                   break;
               case 'c':
                   rptCount = atoi (argv[++i]);
                   if (rptCount < 10) {
                       rptCount = 10;
                   }
                   break;
               default:
                   fprintf (stderr, "Invalid option %s\n", argv[i]);
                   printUsage ();
                   exit (-1);
            }
        }
        ++i;
    }

    list = initList (rptCount);

    matrix = readMatrix (matrixFile);
    if (matrix == NULL) {
        fprintf (stderr, "Error reading matrix\n");
        exit (-1);
    }

    dbLib = openLib (dbFile, swType == WOZNIAK);
    queryLib = openLib (queryFile, 0);

    querySeq = nextSeq (queryLib, &queryLen);
    if (queryLen == 0) {
        fprintf (stderr, "Empty query sequence\n");
        exit (-1);
    }

    printf ("%s vs %s\n", queryFile, dbFile);
    printf ("Matrix: %s, Init: %d, Ext: %d\n\n", 
            matrixFile, options.gapInit, options.gapExt);

    _ftime(&startTime);

    swData = (swFuncs[swType].init) (querySeq, queryLen, matrix);

    (swFuncs[swType].scan) (querySeq, queryLen, dbLib, swData, &options, list);

    (swFuncs[swType].done) (swData);

    _ftime(&endTime);

    printResults (list);

    printf ("\n");
    printf ("%d residues in query string\n", queryLen);
    printf ("%d residues in %d library sequences\n", 
            dbLib->residues, dbLib->sequences);

    startTimeSec = startTime.time + startTime.millitm / 1000.0;
    endTimeSec = endTime.time + endTime.millitm / 1000.0;
    printf ("Scan time: %6.3f (%s implementation)\n", 
            endTimeSec - startTimeSec,
            SW_IMPLEMENATION[swType]);
  
    closeLib (queryLib);
    closeLib (dbLib);

    free (matrix);

    freeList (list);
}
#endif

SCORE_LIST *
initList (int count)
{
    int i;

    SCORE_LIST *hdr;
    SCORE_NODE *list;
    SCORE_NODE *prev;

    hdr = (SCORE_LIST *) malloc (sizeof (SCORE_LIST));
    if (hdr == NULL) {
        fprintf (stderr, "Cannot allocate storage for score header\n");
        exit (-1);
    }

    list = (SCORE_NODE *) calloc (count, sizeof (SCORE_NODE));
    if (list == NULL) {
        fprintf (stderr, "Cannot allocate storage for scores\n");
        exit (-1);
    }

    /* initialize the scores list */
    hdr->minScore = 0;
    hdr->first = NULL;
    hdr->last = NULL;
    hdr->free = list;
    hdr->buffer = list;

    prev = NULL;
    for (i = 0; i < count; ++i) {
        list[i].name[0] = '\0';
        list[i].score = 0;

        if (i == 0) {
            list[i].prev = NULL;
        } else {
            list[i].prev = &list[i-1];
        }

        if (i == count - 1) {
            list[i].next = NULL;
        } else {
            list[i].next = &list[i+1];
        }
    }

    return hdr;
}

void freeList (SCORE_LIST *list)
{
    free (list->buffer);
    free (list);
}

int insertList (SCORE_LIST *list, int score, char *name)
{
    SCORE_NODE *node;
    SCORE_NODE *ptr = list->first;

    if (list->free != NULL) {
        node = list->free;
        list->free = list->free->next;
    } else if (score > list->last->score) {
        node = list->last;
        list->last = node->prev;
        list->last->next = NULL;
    } else {
        /* should never happen */
        return list->minScore + 1;
    }

    strncpy (node->name, name, MAX_SCORE_NAME);
    node->name[MAX_SCORE_NAME - 1] = '\0';
    node->score = score;

    while (ptr && ptr->score >= score) {
        ptr = ptr->next;
    }

    if (list->first == NULL) {
        list->first = node;
        list->last = node;
        node->prev = NULL;
        node->next = NULL;
    } else if (ptr == NULL) {
        node->prev = list->last;
        node->next = NULL;
        node->prev->next  = node;
        list->last = node;
    } else {
        node->prev = ptr->prev;
        node->next = ptr;

        if (node->prev == NULL) {
            list->first = node;
        } else {
            node->prev->next = node;
        }
        ptr->prev = node;
    }

    if (list->free == NULL) {
        list->minScore = list->last->score + 1;
    }

    return list->minScore;
}

void printResults (SCORE_LIST *list)
{
    SCORE_NODE *ptr = list->first;

    printf ("Score  Description\n");

    while (ptr) {
        printf ("%5d  %s\n", ptr->score, ptr->name);
        ptr = ptr->next;
    }
}

