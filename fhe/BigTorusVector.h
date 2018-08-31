#ifndef FHE_BIGTORUSVECTOR_H
#define FHE_BIGTORUSVECTOR_H

#include "BigTorus.h"
#include "commons.h"
#include <cstdint>
#include <memory>

/**
 * A vector of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusVector {
public:
    const BigTorusParams &btp; // a ref to the parameters it was initialized with (only for destructor purposes)
    UINT64 const length;
    UINT64 *const limbs_end;

    NO_COPY(BigTorusVector);

    BigTorusVector(UINT64 length, const BigTorusParams &params);

    BigTorusVector(UINT64 length, BigTorusParams &&params) = delete; //cannot pass a temporary

    ~BigTorusVector();

    BigTorusRef getAT(UINT64 i) const; ///< coef a[i] as a BigTorus
    BigTorusRef getAT(UINT64 i); ///< coef a[i] as a BigTorus
};


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  vector of BigTorus:     length*torus_limbs*UINT64
 */
void serializeBigTorusVectorContent(std::ostream &out, const BigTorusVector &value);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  vector of BigTorus:     length*torus_limbs*UINT64
 */
void deserializeBigTorusVectorContent(std::istream &in, BigTorusVector &reps);

/** serialize:
 *  param:        BigTorusParams
 *  length:       UINT64
 *  value:        BigTorusVectorContent
 */
void serializeBigTorusVector(std::ostream &out, const BigTorusVector &value);

/** serialize:
 *  param:          BigTorusParams
 *  length:         UINT64
 *  content:        BigTorusVectorContent
 */
std::shared_ptr<BigTorusVector> deserializeBigTorusVector(std::istream &in);



/**
 * A matrix of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusMatrix {
public:
    UINT64 *limbs_end; //<end of first bigtorus limb pointer
    UINT64 rows;
    UINT64 cols;
    const BigTorusParams &params;

    BigTorusRef operator()(int i, int j) {
        return BigTorusRef(
                limbs_end + (i * cols + j) * params.torus_limbs,
                params);
    }

    BigTorusRef operator()(int i, int j) const {
        return BigTorusRef(
                limbs_end + (i * cols + j) * params.torus_limbs,
                params);
    };

    BigTorusMatrix(UINT64 rows, UINT64 cols, const BigTorusParams &params);

    ~BigTorusMatrix();
};


void zero(const BigTorusVector &v);

void fixp_add(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits = NA);

void fixp_sub(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits = NA);

/** @brief compute a vector equal to out= oui-a.in */
void subMul128(BigTorusVector &out, __int128 a, const BigTorusVector &in, const UINT64 out_limb_prec);


/** @brief compute a vector equal to out= oui-a.in */
void subMul64(BigTorusVector &out, int64_t a, const BigTorusVector &in, const UINT64 out_limb_prec);


/** @brief copy "in" in "out" */
void copy(BigTorusVector &out, const BigTorusVector &in, UINT64 out_limbs_prec);

/** @brief generate a random bigtorusvector */
void random(BigTorusVector &out, UINT64 out_limbs_prec);

/** @brief add a noise to out  of alpha bits*/
void add_noise(BigTorusVector &out, UINT64 alpha_bits, UINT64 out_limbs_prec);

void sub(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec);

void add(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec);

/** @brief left shift the all torus by exactly shift_bits */
void lshift(BigTorusVector &out, const BigTorusVector &in, int64_t shift_bits);

/** @brief multiply all torus elems by multiplier */
void public_scale(BigTorusVector &out, const BigTorusVector &in, int64_t multiplier);


#endif //FHE_BIGTORUSVECTOR_H

