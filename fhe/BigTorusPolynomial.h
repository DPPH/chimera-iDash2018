#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"

class BigTorusPolynomial : public BigTorusVector {
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};

/** @brief generate a constant polynomial equal to in */
void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec);

#endif //FHE_BIGTORUSPOLYNOMIAL_H
