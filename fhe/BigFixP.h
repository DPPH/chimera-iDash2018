#ifndef FHE_BIGFP_H
#define FHE_BIGFP_H


#include <cstdint>
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>

#include "gmp.h"

const int64_t NA = 1ul << 63ul;
const int64_t BITS_PER_LIMBS = 64;
const int64_t BYTES_PER_LIMBS = BITS_PER_LIMBS / 8;

inline void assert_dramatically(bool condition, const std::string &message = "") {
    if (!condition) {
        std::cerr << "ASSERT_FAILED: " << message << std::endl;
        abort();
    }
}

inline void assert_weakly(bool condition, const std::string &message = "") {
    if (!condition) {
        std::cerr << "ASSERT_WARNING: " << message << std::endl;
    }
}

class BigTorusParams {
public:
    uint64_t torus_limbs; ///< total number of limbs

    BigTorusParams(uint64_t torus_limbs);
};

static_assert(sizeof(BigTorusParams)==8, "Bad compiler");

class BigFixPParams {
public:
    BigTorusParams torus_params; ///< total number of limbs
    int64_t plaintext_expo;      ///< plaintext exponent
    int64_t level_expo;

    BigFixPParams(const BigTorusParams &torus_params, int64_t plaintext_expo, int64_t level_expo);
    ///< level exponent
};

static_assert(sizeof(BigFixPParams)==24, "Bad compiler");

class BigFixP {
public:
    uint64_t *limbs_raw;
    const BigFixPParams * const params;

    BigFixP(const BigFixPParams* params):
            params(params)
    {
        limbs_raw = new uint64_t[params->torus_params.torus_limbs];
    }

    ~BigFixP() {
        delete[] limbs_raw;
    }
};

class BigFixPVector {
public:
    uint64_t length;
    uint64_t *limbs_raw;
    BigFixPParams *params;
};

class BigFixPMatrix {
public:
    uint64_t rows;
    uint64_t cols;
    uint64_t *limbs_raw;
    BigFixPParams *params;
};

class BigTorus {
public:
    uint64_t *limbs_raw;
    const BigTorusParams *const  params;

    BigTorus(const BigTorusParams* params):
    params(params)
    {
        limbs_raw = new uint64_t[params->torus_limbs];
    }

    ~BigTorus() {
        delete[] limbs_raw;
    }
};

class BigTorusVector {
public:
    uint64_t length;
    uint64_t *limbs_raw;
    BigFixPParams *params;
};

class BigTorusMatrix {
public:
    uint64_t rows;
    uint64_t cols;
    uint64_t *limbs_raw;
    BigFixPParams *params;
};

class BigTorusRef {
public:
    uint64_t *limbs_raw;
    const BigTorusParams * const params;

    BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params);
    BigTorusRef(const BigTorus& torus);
    BigTorusRef(BigTorus& torus);
};


class BigFixPRef {
public:
    uint64_t *limbs_raw;
    const BigFixPParams *const params;
    BigFixPRef(const BigFixP& a);
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

void add(BigFixP &reps, const BigFixP &a, const BigFixP &b, int64_t out_precision_bits = NA);

#endif //FHE_BIGFP_H
