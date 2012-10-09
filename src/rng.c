
/* This code taken from the GNU Scientific Library and adapted for use in
 * fastq-tools. It is subject to the following licence.
 */

/* This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA

   Original implementation was copyright (C) 1997 Makoto Matsumoto and
   Takuji Nishimura. Coded by Takuji Nishimura, considering the
   suggestions by Topher Cooper and Marc Rieffel in July-Aug. 1997, "A
   C-program for MT19937: Integer version (1998/4/6)"

   This implementation copyright (C) 1998 Brian Gough. I reorganized
   the code to use the module framework of GSL.  The license on this
   implementation was changed from LGPL to GPL, following paragraph 3
   of the LGPL, version 2.

   Update:

   The seeding procedure has been updated to match the 10/99 release
   of MT19937.

   Update:

   The seeding procedure has been updated again to match the 2002
   release of MT19937

   The original code included the comment: "When you use this, send an
   email to: matumoto@math.keio.ac.jp with an appropriate reference to
   your work".

   Makoto Matsumoto has a web page with more information about the
   generator, http://www.math.keio.ac.jp/~matumoto/emt.html.

   The paper below has details of the algorithm.

   From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
   623-dimensionally equidistributerd uniform pseudorandom number
   generator". ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1 (Jan. 1998), Pages 3-30

   You can obtain the paper directly from Makoto Matsumoto's web page.

   The period of this generator is 2^{19937} - 1.

*/

#include "rng.h"
#include "common.h"
#include <stdlib.h>
#include <assert.h>


/* default seed */
unsigned long int default_seed = 4357;


#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;   

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;   


static const unsigned long RNG_MIN = 0;
static const unsigned long RNG_MAX = 0xffffffffUL;


struct rng_t_
{
    unsigned long mt[N];
    int mti;
};


static inline unsigned long mt_get(rng_t* state)
{
  unsigned long k ;
  unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  if (state->mti >= N)
    {   /* generate N words at one time */
      int kk;

      for (kk = 0; kk < N - M; kk++)
        {
          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
      for (; kk < N - 1; kk++)
        {
          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

      {
        unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
      }

      state->mti = 0;
    }

  /* Tempering */
  
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);

  state->mti++;

  return k;
}

double mt_get_double(rng_t* state)
{
  return mt_get(state) / 4294967296.0;
}

void mt_set(rng_t* state, unsigned long int s)
{
  int i;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  state->mt[0]= s & 0xffffffffUL;

  for (i = 1; i < N; i++)
    {
      /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

      state->mt[i] =
        (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
      
      state->mt[i] &= 0xffffffffUL;
    }

  state->mti = i;
}

#if 0 // these two seeding procedures are not used
static void mt_1999_set(rng_t* state, unsigned long int s)
{
  int i;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  /* This is the October 1999 version of the seeding procedure. It
     was updated by the original developers to avoid the periodicity
     in the simple congruence originally used.

     Note that an ANSI-C unsigned long integer arithmetic is
     automatically modulo 2^32 (or a higher power of two), so we can
     safely ignore overflow. */

#define LCG(x) ((69069 * x) + 1) &0xffffffffUL

  for (i = 0; i < N; i++)
    {
      state->mt[i] = s & 0xffff0000UL;
      s = LCG(s);
      state->mt[i] |= (s &0xffff0000UL) >> 16;
      s = LCG(s);
    }

  state->mti = i;
}

/* This is the original version of the seeding procedure, no longer
   used but available for compatibility with the original MT19937. */

static void mt_1998_set(rng_t* state, unsigned long int s)
{
  int i;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  state->mt[0] = s & 0xffffffffUL;

#define LCG1998(n) ((69069 * n) & 0xffffffffUL)

  for (i = 1; i < N; i++)
    state->mt[i] = LCG1998 (state->mt[i - 1]);

  state->mti = i;
}
#endif

rng_t* fastq_rng_alloc()
{
    rng_t* rng = malloc_or_die(sizeof(rng_t));
    mt_set(rng, default_seed);
    return rng;
}


void fastq_rng_free(rng_t* rng)
{
    free(rng);
}

void fastq_rng_seed(rng_t* rng, unsigned long seed)
{
    mt_set(rng, seed);
}

unsigned long fastq_rng_uniform_int(rng_t* rng, unsigned long k)
{
    unsigned long scale = (RNG_MAX - RNG_MIN) / k;
    assert(scale > 0);
    unsigned long r;

    do {
        r = (mt_get(rng) - RNG_MIN) / scale;
    } while (r >= k);

    return r;
}

