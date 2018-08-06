#ifndef FHE_BIGTORUSVECTOR_H
#define FHE_BIGTORUSVECTOR_H


#include "BigFixP.h"

/**
 * A vector of BigTorus (with shared params), owns its limbs memory space
 */
class BigTorusVector {
protected:
    const BigTorusParams &btp; // a ref to the parameters it was initialized with (only for destructor purposes)
public:
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
    uint64_t *limbs_raw;
    uint64_t rows;
    uint64_t cols;
    BigTorusParams *params;
};


#endif //FHE_BIGTORUSVECTOR_H
