
#include "hash.h"
#include "common.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static const size_t INITIAL_TABLE_SIZE = 128;
static const double MAX_LOAD = 0.75;


/*
 * Paul Hsieh's SuperFastHash
 * http://www.azillionmonkeys.com/qed/hash.html
 */


#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
    || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
        +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

static uint32_t hash(const char * data, size_t len) {
    uint32_t hash = len, tmp;
    int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}



static void rehash(hash_table* T, size_t new_n);
static void clear_hash_table(hash_table*);



hash_table* create_hash_table()
{
    hash_table* T = malloc_or_die(sizeof(hash_table));
    T->A = malloc_or_die(INITIAL_TABLE_SIZE * sizeof(hashed_value*));
    memset(T->A, 0, INITIAL_TABLE_SIZE * sizeof(hashed_value*));
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;

    return T;
}


void destroy_hash_table(hash_table* T)
{
    if (T != NULL) {
        clear_hash_table(T);
        free(T->A);
        free(T);
    }
}



void clear_hash_table(hash_table* T)
{
    hashed_value* u;
    size_t i;
    for (i = 0; i < T->n; i++) {
        while (T->A[i]) {
            u = T->A[i]->next;
            free(T->A[i]->value);
            free(T->A[i]);
            T->A[i] = u;
        }
    }

    T->m = 0;
}


static void insert_without_copy(hash_table* T, hashed_value* V)
{
    uint32_t h = hash(V->value, V->len) % T->n;
    V->next = T->A[h];
    T->A[h] = V;
    T->m++;
}



static void rehash(hash_table* T, size_t new_n)
{
    hash_table U;
    U.n = new_n;
    U.m = 0;
    U.max_m = U.n * MAX_LOAD;
    U.A = malloc_or_die(U.n * sizeof(hashed_value*));
    memset(U.A, 0, U.n * sizeof(hashed_value*));

    hashed_value *j, *k;
    size_t i;
    for (i = 0; i < T->n; i++) {
        j = T->A[i];
        while (j) {
            k = j->next;
            insert_without_copy(&U, j);
            j = k;
        }
        T->A[i] = NULL;
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = U.max_m;
}


void inc_hash_table(hash_table* T, const char* value, size_t len)
{
    if (T->m >= T->max_m) rehash(T, T->n * 2);

    uint32_t h = hash(value, len) % T->n;

    hashed_value* u = T->A[h];

    while (u) {
        if (u->len == len && memcmp(u->value, value, len) == 0) {
            u->count++;
            return;
        }

        u = u->next;
    }

    u = malloc_or_die(sizeof(hashed_value));
    u->value = malloc_or_die(len);
    memcpy(u->value, value, len);
    u->len = len;
    u->count = 1;

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}



hashed_value** dump_hash_table(hash_table* T)
{
    hashed_value** D = malloc_or_die(T->m * sizeof(hashed_value*));

    hashed_value* u;
    size_t i, j;
    for (i = 0, j = 0; i < T->n; i++) {
        u = T->A[i];
        while (u) {
            D[j++] = u; 
            u = u->next;
        }
    }

    return D;
}



