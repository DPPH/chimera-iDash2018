#ifndef FHE_BIGFIXPVECTOR_H
#define FHE_BIGFIXPVECTOR_H


#include "BigFixP.h"
#include "BigTorusVector.h"


class BigFixPVector : public BigTorusVector {
public:
    const BigFixPParams &bfp; ///< fixed point params of the vector (cast of btp)
    BigFixPRef getAF(uint64_t i); ///< get i-th elem as fixed point
    BigFixPRef getAF(uint64_t i) const; ///< get i-th elem as fixed point

    BigFixPVector(uint64_t length, const BigFixPParams &params);

    ~BigFixPVector();
};

class BigFixPMatrix {
public:
    uint64_t *limbs_raw;
    const uint64_t rows;
    const uint64_t cols;
    const BigFixPParams *const params;

    BigFixPRef operator()(uint64_t i, uint64_t j);

    BigFixPMatrix(uint64_t rows, uint64_t cols, const BigFixPParams *params);

    ~BigFixPMatrix();
};

void clear(BigFixPVector &reps);

void add(BigFixPVector &reps, const BigFixPVector &a, const BigFixPVector &b, uint64_t out_precision_bits = NA);

void sub(BigFixPVector &reps, const BigFixPVector &a, const BigFixPVector &b, uint64_t out_precision_bits = NA);

#endif //FHE_BIGFIXPVECTOR_H
