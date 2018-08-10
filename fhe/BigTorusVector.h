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



#endif //FHE_BIGTORUSVECTOR_H
