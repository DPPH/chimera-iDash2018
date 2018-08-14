#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"
#include "BigReal.h"
#include "BigComplex.h"

class BigTorusPolynomial : public BigTorusVector {
public:
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};

/** @brief generate a constant polynomial equal to in */
void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec);

/** @brief int-polynomial torus-polynomial external product  */
void naive_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, UINT64 out_limb_prec);


/** @brief int-polynomial torus-polynomial external product using fft  */
void fft_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, const UINT64 bits_a,
                          UINT64 out_limb_prec);

/** @brief polynomial-polynomial torus-polynomial internal product using fft  */
void fft_internal_product(BigTorusPolynomial &out, const BigTorusPolynomial &a, const BigTorusPolynomial &b,
                          UINT64 out_limb_prec);

/** @brief int-polynomial torus-polynomial external product using fft  (TODO: a voir comment definir a)*/
void fft_semi_external_product(BigTorusPolynomial &out, const BigComplex *ca, const BigTorusPolynomial &b,
                               const UINT64 bits_a, UINT64 out_limb_prec);


#endif //FHE_BIGTORUSPOLYNOMIAL_H
