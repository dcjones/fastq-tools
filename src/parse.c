/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include "parse.h"
#include "common.h"
#include <stdlib.h>
#include <ctype.h>


static const size_t init_str_size  = 128;
static const size_t fastq_buf_size = 4096;

static void fastq_alloc_str(str_t* s)
{
    s->s = malloc_or_die(init_str_size);
    s->s[0] = '\0';
    s->n = 0;
    s->size = init_str_size;
}


static void fastq_expand_str(str_t* s)
{
    s->size *= 2;
    realloc_or_die(s->s, s->size);
}


seq_t* fastq_alloc_seq()
{
    seq_t* seq = malloc_or_die(sizeof(seq_t));
    fastq_alloc_str(&seq->id1);
    fastq_alloc_str(&seq->seq);
    fastq_alloc_str(&seq->id2);
    fastq_alloc_str(&seq->qual);

    return seq;
}


void fastq_free_seq(seq_t* seq)
{
    free(seq->id1.s);
    free(seq->seq.s);
    free(seq->id2.s);
    free(seq->qual.s);
    free(seq);
}


typedef enum
{
    STATE_EOF,
    STATE_ERROR,
    STATE_ID1,
    STATE_SEQ,
    STATE_ID2,
    STATE_QUAL

} fastq_state;


fastq_t* fastq_open(FILE* f)
{
    fastq_t* fqf = malloc_or_die(sizeof(fastq_t));
    or_die((int)((fqf->file = gzdopen(fileno(f), "rb")) != NULL),
           "Can not open gzip file.");
    
    fqf->state = STATE_ID1;
    fqf->buf = malloc_or_die(fastq_buf_size);
    fqf->buf[0] = '\0';
    fqf->c = fqf->buf;

    return fqf;
}


void fastq_close(fastq_t* fqf)
{
    gzclose(fqf->file);
    free(fqf->buf);
    free(fqf);
}


void fastq_refill(fastq_t* f)
{
    int errnum;
    const char* errmsg;

    int n = gzread(f->file, f->buf, fastq_buf_size - 1);

    if (n <= 0) {
        if (gzeof(f->file)) {
            f->state = STATE_EOF;
            n = 0;
        }
        else {
            errmsg = gzerror(f->file, &errnum);
            fprintf(stderr, "I/O error: %s\n", errmsg);
            exit(1);
        }
    }

    f->buf[n] = '\0';
    f->c = f->buf;
}


void fastq_get_line(fastq_t* f, str_t* s)
{
    size_t i = 0;

    if (f->state == STATE_EOF) goto fastq_get_line_done;

    while (1) {
        switch (*f->c) {
            case '\0':
                fastq_refill(f);
                if (f->state == STATE_EOF) goto fastq_get_line_done;
                break;

            case '\r':
                f->c++;
                break;

            case '\n':
                goto fastq_get_line_done;

            default:
                while (s->size < i + 2) {
                    fastq_expand_str(s);
                }
                if (s) s->s[i++] = *f->c;
                f->c++;
        }

    }

fastq_get_line_done:
    if (s) {
        s->s[i] = '\0';
        s->n = i;
    }
}



int fastq_next(fastq_t* f, seq_t* seq)
{
    if (f->state == STATE_EOF) return 0;

    while (1) {

        /* read more, if needed */
        if (*f->c == '\0' ) {
            fastq_refill(f);
            if (f->state == STATE_EOF) return 0;
            continue;
        }

        /* skip over leading whitespace */
        else if (isspace(*f->c)) {
            /* do nothing */
        }

        /* skip comments */
        /*
        else if (*f->c == ';') {
            fastq_get_line(f, NULL);
            if (f->state == STATE_EOF) return 0;
        }
        */

        /* read id1 */
        else if (f->state == STATE_ID1) {
            if (*f->c == '@' || *f->c == '>') {
                f->c++;
                fastq_get_line(f, &seq->id1);
                if (f->state == STATE_EOF) return 0;

                f->state = STATE_SEQ;
            }
            else {
                fprintf(stderr,
                        "Malformed FASTQ file: expecting an '@' or '>', saw a '%c'\n",
                        *f->c);
                exit(1);
            }
        }

        /* read sequence */
        else if (f->state == STATE_SEQ) {
            fastq_get_line(f, &seq->seq);
            if (f->state == STATE_EOF) return 0;

            f->state = STATE_ID2;
        }

        /* read id2 */
        else if (f->state == STATE_ID2) {
            if (*f->c == '+') {
                f->c++;
                fastq_get_line(f, &seq->id2);
                if (f->state == STATE_EOF) return 0;

                f->state = STATE_QUAL;
            }
            else {
                /* fasta style entry */
                seq->id2.s[0]  = '\0';
                seq->qual.s[0] = '\0';

                f->state = STATE_ID1;
                break;
            }
        }

        /* read quality string */
        else if (f->state == STATE_QUAL) {
            fastq_get_line(f, &seq->qual);
            if (f->state == STATE_EOF) return 1;

            f->state = STATE_ID1;
            break;
        }

        else {
            fputs("Inexplicable error in fastq parser.\n", stderr);
            exit(1);
        }

        f->c++;
    }

    return 1;
}


void fastq_print(FILE* fout, seq_t* seq)
{
    /* FASTQ */
    if (seq->qual.n > 0) {
        fprintf(fout, "@%s\n%s\n+%s\n%s\n",
                      seq->id1.s,
                      seq->seq.s,
                      seq->id2.s,
                      seq->qual.s );
    }

    /* FASTA */
    else {
        fprintf(fout, ">%s\n%s\n",
                      seq->id1.s,
                      seq->seq.s );
    }
}


