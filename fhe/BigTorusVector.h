#ifndef FHE_BIGTORUSVECTOR_H
#define FHE_BIGTORUSVECTOR_H

#include "BigTorus.h"

/**
 * A vector of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusVector {
public:
    const BigTorusParams &btp; // a ref to the parameters it was initialized with (only for destructor purposes)
    uint64_t const length;
    uint64_t *const limbs;

    NO_COPY(BigTorusVector);

    BigTorusVector(uint64_t length, const BigTorusParams &params);

    ~BigTorusVector();

    BigTorusRef getAT(uint64_t i) const; ///< coef a[i] as a BigTorus
    BigTorusRef getAT(uint64_t i); ///< coef a[i] as a BigTorus
};

/**
 * A matrix of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusMatrix {
public:
    uint64_t *limbs;
    uint64_t rows;
    uint64_t cols;
    const BigTorusParams *params;

    BigTorusRef operator()(int i, int j) { return BigTorusRef(limbs + (i * cols + j) * params->torus_limbs, params); }

    BigTorusRef operator()(int i, int j) const {
        return BigTorusRef(limbs + (i * cols + j) * params->torus_limbs, params);
    };

    BigTorusMatrix(uint64_t rows, uint64_t cols, const BigTorusParams *params);

    ~BigTorusMatrix();
};

void zero(const BigTorusVector &v);

void fixp_add(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, uint64_t out_precision_bits = NA);

void fixp_sub(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, uint64_t out_precision_bits = NA);



#endif //FHE_BIGTORUSVECTOR_H
