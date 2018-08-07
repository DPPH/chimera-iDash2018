#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"

class BigTorusPolynomial : public BigTorusVector {
    BigTorusPolynomial(uint64_t N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};


#endif //FHE_BIGTORUSPOLYNOMIAL_H
