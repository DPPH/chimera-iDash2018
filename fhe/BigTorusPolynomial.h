#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"

class BigTorusPolynomial : public BigTorusVector {
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};


#endif //FHE_BIGTORUSPOLYNOMIAL_H
