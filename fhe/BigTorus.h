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
    const UINT64 torus_limbs; ///< total number of limbs
    const int64_t plaintext_expo;  ///< plaintext exponent
    const int64_t level_expo;      ///< level exponent

    BigTorusParams(UINT64 torus_limbs, int64_t plaintext_expo = 0, int64_t level_expo = 0);
};

static_assert(sizeof(BigTorusParams) == 24, "Bad compiler");

/**
 * A standalone BigTorus (owns its limbs memory space)
 */
class BigTorus {
public:
    UINT64 *limbs_end;            ///< limbs end ptr
    const BigTorusParams &params; ///< parameters

    /** constructs new BigTorus */
    BigTorus(const BigTorusParams &params) :
            params(params) {
        UINT64 *limbs = new UINT64[params.torus_limbs];
        memset(limbs, 0, params.torus_limbs * sizeof(UINT64));
        limbs_end = limbs + params.torus_limbs;

    }

    BigTorus(BigTorusParams &&params) = delete; //cannot pass a temporary


    /** deletes a bigtorus */
    ~BigTorus() {
        delete[] (limbs_end - params.torus_limbs);
    }

    BigTorus(const BigTorus &) = delete;

    void operator=(const BigTorus &)= delete;
};

/**
 * A weak reference to a BigTorus (does not own memory)
 */
class BigTorusRef {
public:
    UINT64 *const limbs_end;
    const BigTorusParams &params;

    BigTorusRef(UINT64 *limbs_end, const BigTorusParams &params);

    BigTorusRef(const BigTorus &torus);

    BigTorusRef(BigTorus &torus);

    BigTorusRef(BigTorus &&torus) = delete; //cannot pass a temporary
};

/** multiply a bigTorus by an integer (modulo 1) */
void bigTorusRawScaleU64(UINT64 *limbs_end, int64_t coef, UINT64 nblimbs);

/** multiply a bigTorus by an integer (modulo 1) */
void bigTorusScaleU64(const BigTorusRef &x, int64_t coef);

/** copy the limb_precision most significant bits from x to dest */
void copy(BigTorusRef dest, const BigTorusRef &x, UINT64 limb_precision = NA);

/** make dest zero */
void zero(BigTorusRef dest);

/** make the limb_precision most significant bits of dest uniformly random */
void random(BigTorusRef dest, UINT64 limb_precision = NA);

/** add a noise of magnitude 2^-alpha (cutoff distribution for now) */
void add_noise(BigTorusRef dest, UINT64 alpha_bits, UINT64 limb_precision = NA);

struct bigtorus_addsub_params {
    UINT64 asize;
    UINT64 bsize;
    UINT64 dsize;
    UINT64 limb_precision;
};

/** prepare an addition or a subtraction (with the provided precision) */
bigtorus_addsub_params prepare_addsub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b,
                                      UINT64 limb_precision = NA);

/** add two bigtorus (prepared) */
void add_prep(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, const bigtorus_addsub_params &params);

/** add two bigtorus (prepared) */
void sub_prep(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, const bigtorus_addsub_params &params);

/** add two bigtorus (with the provided precision, wrapper function) */
void add(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, UINT64 limb_precision = NA);

/** add two bigtorus (with the provided precision, wrapper function) */
void sub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, UINT64 limb_precision = NA);

/** @brief compute out=in*/
void neg(BigTorusRef out, const BigTorusRef &in, UINT64 limb_precision = NA);

/** @brief compute the bigtorus equal to out= out-a.in */
void subMulS128(BigTorusRef out, __int128 a, const BigTorusRef &in, UINT64 out_limb_prec);

/** @brief compute the bigtorus equal to out= out-a.in */
void subMulS64(BigTorusRef out, int64_t a, const BigTorusRef &in, UINT64 out_limb_prec);

/** @brief compute out = a * in */
void mulS64(BigTorusRef out, int64_t a, const BigTorusRef &in, UINT64 limb_prec);

void to_torus(BigTorusRef reps, const NTL::RR &a);

double log2Diff(const BigTorusRef &a, const BigTorusRef &b);

NTL::RR to_RR(const BigTorusRef &a);

// ***********************

struct fixp_add_params {
    int64_t ia_limbs;   //total limbs in a
    int64_t ib_limbs;   //total limbs in b
    int64_t ic_limbs;   //total limbs in c
    int64_t oc_limbs;   //significant limbs in c's output
    int64_t sa;         // shift in bits to apply to a
    int64_t sb;         // shift in bits to apply to b
    UINT64 *tmp_a; //buffer of size oc_limbs+1
    UINT64 *tmp_b; //buffer of size oc_limbs+1

    fixp_add_params();
};


void fixp_prepareAdd(fixp_add_params &out,
                     const BigTorusParams &pres, const BigTorusParams &pa, const BigTorusParams &pb,
                     int64_t out_precision_bits
);

void fixp_releaseAdd(fixp_add_params &out);

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
void fixp_special_add_lshift(UINT64 *out, UINT64 *a,
                             int64_t out_limbs, int64_t a_limbs, int64_t bit_shift);

/**
 * Performs the raw fixed point addition on the torus buffers a and b
 *
 * @param reps torus buffer for the result
 * @param a torus buffer of a
 * @param b torus buffer of b
 * @param params the addition parameters
 */
void fixp_raw_add(UINT64 *reps, UINT64 *a, UINT64 *b, const fixp_add_params &params);

/**
 * Performs the raw fixed point subtraction on the torus buffers a and b
 *
 * @param reps torus buffer for the result
 * @param a torus buffer of a
 * @param b torus buffer of b
 * @param params the addition parameters
 */
void fixp_raw_sub(UINT64 *reps, UINT64 *a, UINT64 *b, const fixp_add_params &params);

/**
 * zero
 */
void fixPRawClear(UINT64 *reps, const UINT64 limbs_size);

void fixp_add(BigTorusRef reps, const BigTorusRef &a, const BigTorusRef &b, UINT64 out_precision_bits = NA);

void fixp_sub(BigTorusRef reps, const BigTorusRef &a, const BigTorusRef &b, UINT64 out_precision_bits = NA);


#endif //FHE_BIGTORUS_H
