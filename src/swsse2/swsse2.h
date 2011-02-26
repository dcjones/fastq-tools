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

#ifndef INCLUDE_SWSSE2_H
#define INCLUDE_SWSSE2_H

typedef void SW_DATA;

#define ALPHA_SIZE 23

extern const char AMINO_ACIDS[ALPHA_SIZE];
extern const int AMINO_ACID_VALUE[256];

#define SHORT_BIAS 32768

typedef struct {
    short gapInit;
    short gapExt;
    int   threshold;
} SEARCH_OPTIONS;

#define MAX_SCORE_NAME 64
typedef struct SCORE_STRUCT {
    int score;
    char name[MAX_SCORE_NAME];
    struct SCORE_STRUCT *prev;
    struct SCORE_STRUCT *next;
} SCORE_NODE;

typedef struct {
    int minScore;
    SCORE_NODE *first;
    SCORE_NODE *last;
    SCORE_NODE *free;
    void *buffer;
} SCORE_LIST;

int insertList (SCORE_LIST *list, int score, char *name);

#endif /* INCLUDE_SWSSE2_H */
