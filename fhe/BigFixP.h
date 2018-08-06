#ifndef FHE_BIGFP_H
#define FHE_BIGFP_H


#include <cstdint>
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <gmp.h>

#include "commons.h"
#include "BigTorus.h"


class BigFixPParams : public BigTorusParams {
public:
    int64_t plaintext_expo;      ///< plaintext exponent
    int64_t level_expo; ///< level exponent

    BigFixPParams(const BigTorusParams &torus_params, int64_t plaintext_expo, int64_t level_expo);
};

static_assert(sizeof(BigFixPParams) == 24, "Bad compiler");

class BigFixP {
public:
    uint64_t *limbs_raw;
    const BigFixPParams *const params;

    BigFixP(const BigFixPParams *params) :
            params(params) {
        limbs_raw = new uint64_t[params->torus_limbs];
    }

    ~BigFixP() {
        delete[] limbs_raw;
    }
};




class BigFixPRef {
public:
    uint64_t *limbs_raw;
    const BigFixPParams *const params;

    BigFixPRef(const BigFixP &a);

    BigFixPRef(uint64_t *limbs, const BigFixPParams *params);
};


// ***********************

class BigFixPAddParams {
public:
    int64_t ia_limbs;   //total limbs in a
    int64_t ib_limbs;   //total limbs in b
    int64_t ic_limbs;   //total limbs in c
    int64_t oc_limbs;   //significant limbs in c's output
    int64_t sa;         // shift in bits to apply to a
    int64_t sb;         // shift in bits to apply to b
    uint64_t *tmp_a; //buffer of size oc_limbs+1
    uint64_t *tmp_b; //buffer of size oc_limbs+1

    BigFixPAddParams();
};


void prepareAdd(BigFixPAddParams &out,
                const BigFixPParams &pres, const BigFixPParams &pa, const BigFixPParams &pb,
                int64_t out_precision_bits
);

void releaseAdd(BigFixPAddParams &out);

/**
 * left shift a by bit_shift bits, and write the out_limbs msb of the result to out.
 * Warning: out must have physical size at least out_limbs+1
 *
 * @param out output torus fractional buffer of size out_limbs (physical out_limbs+1)
 * @param a torus fractional buffer of size a_limbs
 * @param out_limbs
 * @param a_limbs
 * @param bit_shift left shift in bits
 */
void special_add_lshift(uint64_t *out, uint64_t *a,
                        int64_t out_limbs, int64_t a_limbs, int64_t bit_shift);

/**
 * Performs the raw fixed point addition on the torus buffers a and b
 *
 * @param reps torus buffer for the result
 * @param a torus buffer of a
 * @param b torus buffer of b
 * @param params the addition parameters
 */
void fixPRawAdd(uint64_t *reps, uint64_t *a, uint64_t *b, const BigFixPAddParams &params);

/**
 * Performs the raw fixed point subtraction on the torus buffers a and b
 *
 * @param reps torus buffer for the result
 * @param a torus buffer of a
 * @param b torus buffer of b
 * @param params the addition parameters
 */
void fixPRawSub(uint64_t *reps, uint64_t *a, uint64_t *b, const BigFixPAddParams &params);

/**
 * zero
 */
void fixPRawClear(uint64_t *reps, const uint64_t limbs_size);

void add(BigFixP &reps, const BigFixP &a, const BigFixP &b, uint64_t out_precision_bits = NA);

void sub(BigFixP &reps, const BigFixP &a, const BigFixP &b, uint64_t out_precision_bits = NA);

void clear(BigFixP &reps);

#endif //FHE_BIGFP_H
