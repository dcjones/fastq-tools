
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "parse.h"
#include "common.h"


static void str_init(str_t* str)
{
    str->n = 0;
    str->size = 128;
    str->s = malloc_or_die(str->size);
    str->s[0] = '\0';
}


static void str_free(str_t* str)
{
    free(str->s);
}


/* Reserve space for `size` more characters. */
static void str_reserve_extra(str_t* str, size_t size)
{
    if (str->n + size > str->size) {
        if (str->n + size > 2 * str->size) {
            str->size = str->n + size;
        }
        else str->size *= 2;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}


/* Copy n characters from c to the end of str. */
static void str_append(str_t* str, char* c, size_t n)
{
    str_reserve_extra(str, n);
    memcpy(str->s + str->n, c, n);
    str->n += n;
    str->s[str->n] = '\0';
}


seq_t* seq_create()
{
    seq_t* seq = malloc_or_die(sizeof(seq_t));
    str_init(&seq->id1);
    str_init(&seq->seq);
    str_init(&seq->id2);
    str_init(&seq->qual);
    return seq;
}


void seq_free(seq_t* seq)
{
    str_free(&seq->id1);
    str_free(&seq->seq);
    str_free(&seq->id2);
    str_free(&seq->qual);
    free(seq);
}


/* This is MurmurHash3. The original C++ code was placed in the public domain
 * by its author, Austin Appleby. */

static inline uint32_t fmix(uint32_t h)
{
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    return h;
}


static inline uint32_t rotl32(uint32_t x, int8_t r)
{
    return (x << r) | (x >> (32 - r));
}


uint32_t murmurhash3(uint32_t seed, const uint8_t* data, size_t len_)
{
    const int len = (int) len_;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    uint32_t c1 = 0xcc9e2d51;
    uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t * blocks = (const uint32_t*) (data + nblocks * 4);

    int i;
    for(i = -nblocks; i; i++)
    {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;

        h1 ^= k1;
        h1 = rotl32(h1, 13);
        h1 = h1*5+0xe6546b64;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

    uint32_t k1 = 0;

    switch(len & 3)
    {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
              k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
    }

    //----------
    // finalization

    h1 ^= len;

    h1 = fmix(h1);

    return h1;
}


static uint32_t seq_hash_seed = 0xc062fb4a;


void seq_hash_set_seed(uint32_t seed)
{
    seq_hash_seed = seed;
}


uint32_t seq_hash(const seq_t* seq)
{
    uint32_t h = seq_hash_seed;
    h = murmurhash3(h, (uint8_t*) seq->id1.s, seq->id1.n);
    h = murmurhash3(h, (uint8_t*) seq->seq.s, seq->seq.n);
    h = murmurhash3(h, (uint8_t*) seq->id2.s, seq->id2.n);
    h = murmurhash3(h, (uint8_t*) seq->qual.s, seq->qual.n);
    return h;
}


static const size_t parser_buf_size = 1000000;


struct fastq_t_
{
    FILE* file;
    size_t readlen;
    char* buf;
    char* next;
    bool linestart;
};



fastq_t* fastq_create(FILE* file)
{
    fastq_t* f = malloc_or_die(sizeof(fastq_t));
    f->file = file;
    f->next = f->buf = malloc_or_die(parser_buf_size);
    f->readlen = 0;
    f->linestart = true;
    return f;
}


void fastq_free(fastq_t* f)
{
    free(f->buf);
    free(f);
}


typedef enum {
    FASTQ_STATE_ID1,  /* Reading ID1. */
    FASTQ_STATE_SEQ,  /* Reading the sequence. */
    FASTQ_STATE_ID2,  /* Reading ID2. */
    FASTQ_STATE_QUAL, /* Reading quality scores. */
} fastq_parser_state_t;


bool fastq_read(fastq_t* f, seq_t* seq)
{
    seq->id1.n = seq->seq.n = seq->id2.n = seq->qual.n = 0;
    fastq_parser_state_t state = FASTQ_STATE_ID1;
    char* end = f->buf + f->readlen;
    do {
        while (f->next < end) {
            /* Consume pointless special characters prefixing IDs */
            if ((state == FASTQ_STATE_ID1 && f->linestart && f->next[0] == '@') ||
                (state == FASTQ_STATE_ID2 && f->linestart && f->next[0] == '+')) {
                f->linestart = false;
                ++f->next;
                continue;
            }

            char* u = memchr(f->next, '\n', end - f->next);
            if (u == NULL) {
                f->linestart = false;
                u = end;
            }
            else f->linestart = true;

            switch (state) {
                case FASTQ_STATE_ID1:
                    str_append(&seq->id1, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_SEQ;
                    break;

                case FASTQ_STATE_SEQ:
                    str_append(&seq->seq, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_ID2;
                    break;

                case FASTQ_STATE_ID2:
                    str_append(&seq->id2, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_QUAL;
                    break;

                case FASTQ_STATE_QUAL:
                    str_append(&seq->qual, f->next, u - f->next);
                    if (f->linestart) {
                        f->next = u + 1;
                        return true;
                    }
                    break;
            }

            f->next = u + 1;
        }

        /* Try to read more. */
        f->readlen = fread(f->buf, 1, parser_buf_size, f->file);
        f->next = f->buf;
        end = f->buf + f->readlen;
    } while (f->readlen);

    return false;
}


void fastq_rewind(fastq_t* f)
{
    rewind(f->file);
    f->next = f->buf;
    f->readlen = 0;
}


void fastq_print(FILE* fout, const seq_t* seq)
{
    fprintf(fout, "@%s\n%s\n+%s\n%s\n",
                  seq->id1.s,
                  seq->seq.s,
                  seq->id2.s,
                  seq->qual.s );
}


