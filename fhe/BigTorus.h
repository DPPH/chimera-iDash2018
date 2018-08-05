#ifndef FHE_BIGTORUS_H
#define FHE_BIGTORUS_H

#include <cstdint>
#include "commons.h"

/**
 * A bigTorus holds a real number modulo 1, with maximal precision N*64 bits.
 * Internally, a BigTorus is an array of N limbs, the least significant one first.
 *
 * The unsigned value of the BigTorus is sum_{i=0}^{N-1} limb_i 2^{-64(N-i)} mod 1.
 *
 * If we need the value of the BigTorus in [-0.5,0.5[, the sign of the BigTorus is
 * therefore the sign of its msb limb.
 */

/**
 * BigTorus parameters
 */
class BigTorusParams {
public:
    uint64_t torus_limbs; ///< total number of limbs

    BigTorusParams(uint64_t torus_limbs);
};

static_assert(sizeof(BigTorusParams) == 8, "Bad compiler");

/**
 * A standalone BigTorus (owns its limbs memory space)
 */
class BigTorus {
public:
    uint64_t *limbs_raw; ///< limbs
    const BigTorusParams *const params; ///< parameters

    /** constructs new BigTorus */
    BigTorus(const BigTorusParams *params) :
            params(params) {
        limbs_raw = new uint64_t[params->torus_limbs];
    }

    /** deletes a bigtorus */
    ~BigTorus() {
        delete[] limbs_raw;
    }

    BigTorus(const BigTorus &) = delete;

    void operator=(const BigTorus &)= delete;
};

/**
 * A vector of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusVector {
public:
    uint64_t *limbs_raw;
    uint64_t length;
    BigTorusParams *params;
};

/**
 * A matrix of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusMatrix {
public:
    uint64_t *limbs_raw;
    uint64_t rows;
    uint64_t cols;
    BigTorusParams *params;
};

/**
 * A weak reference to a BigTorus (does not own memory)
 */
class BigTorusRef {
public:
    uint64_t *limbs_raw;
    const BigTorusParams *const params;

    BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params);

    BigTorusRef(const BigTorus &torus);

    BigTorusRef(BigTorus &torus);
};

/** multiply a bigTorus by an integer (modulo 1) */
void bigTorusRawScale(uint64_t *limbs, int64_t coef, uint64_t nblimbs);

/** multiply a bigTorus by an integer (modulo 1) */
void bigTorusScale(const BigTorusRef &x, int64_t coef);

/** copy the limb_precision most significant bits from x to dest */
void copy(BigTorusRef dest, const BigTorusRef &x, uint64_t limb_precision = NA);

/** copy the limb_precision most significant bits from x to dest */
void random(BigTorusRef dest, uint64_t limb_precision = NA);

/** add two bigtorus (with the provided precision) */
void add(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision = NA);

/** add two bigtorus (with the provided precision) */
void sub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision = NA);

/** add a noise of magnitude 2^-alpha (cutoff distribution for now) */
void add_noise(BigTorusRef dest, uint64_t alpha_bits, uint64_t limb_precision = NA);

void to_torus(BigTorusRef reps, const NTL::RR &a);

NTL::RR to_RR(const BigTorusRef &a);


#endif //FHE_BIGTORUS_H
