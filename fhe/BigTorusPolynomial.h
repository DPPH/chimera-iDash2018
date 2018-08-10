#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"

class BigTorusPolynomial : public BigTorusVector {
public:
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};

/** @brief generate a constant polynomial equal to in */
void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec);

/** @brief int-polynomial torus-polynomial external product  */
void naive_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, UINT64 out_limb_prec);

#endif //FHE_BIGTORUSPOLYNOMIAL_H
