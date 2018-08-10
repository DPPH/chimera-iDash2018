#ifndef FHE_BIGTORUSVECTOR_H
#define FHE_BIGTORUSVECTOR_H

#include "BigTorus.h"

/**
 * A vector of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusVector {
public:
    const BigTorusParams &btp; // a ref to the parameters it was initialized with (only for destructor purposes)
    UINT64 const length;
    UINT64 *const limbs;

    NO_COPY(BigTorusVector);

    BigTorusVector(UINT64 length, const BigTorusParams &params);

    ~BigTorusVector();

    BigTorusRef getAT(UINT64 i) const; ///< coef a[i] as a BigTorus
    BigTorusRef getAT(UINT64 i); ///< coef a[i] as a BigTorus
};

/**
 * A matrix of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusMatrix {
public:
    UINT64 *limbs;
    UINT64 rows;
    UINT64 cols;
    const BigTorusParams *params;

    BigTorusRef operator()(int i, int j) { return BigTorusRef(limbs + (i * cols + j) * params->torus_limbs, params); }

    BigTorusRef operator()(int i, int j) const {
        return BigTorusRef(limbs + (i * cols + j) * params->torus_limbs, params);
    };

    BigTorusMatrix(UINT64 rows, UINT64 cols, const BigTorusParams *params);

    ~BigTorusMatrix();
};

void zero(const BigTorusVector &v);

void fixp_add(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits = NA);

void fixp_sub(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits = NA);

/** @brief compute a vector equal to out= oui-a.in */
void subMul(BigTorusVector &out, __int128 a, const BigTorusVector &in, const UINT64 out_limb_prec);

/** @brief copy "in" in "out" */
void copy(BigTorusVector &out, const BigTorusVector &in, UINT64 out_limbs_prec);

/** @brief generate a random bigtorusvector */
void random(BigTorusVector &out, UINT64 out_limbs_prec);

/** @brief add a noise to out  of alpha bits*/
void add_noise(BigTorusVector &out, UINT64 alpha_bits, UINT64 out_limbs_prec);

void sub(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec);

void add(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec);

#endif //FHE_BIGTORUSVECTOR_H

