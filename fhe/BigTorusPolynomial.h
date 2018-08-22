#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"
#include "BigReal.h"
#include "BigComplex.h"

class BigTorusPolynomial : public BigTorusVector {
public:
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    BigTorusPolynomial(UINT64 N, BigTorusParams &&params) = delete;

    ~BigTorusPolynomial();
};

/** @brief generate a constant polynomial equal to in */
void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec);

/** @brief addition */
void add(BigTorusPolynomial &reps, const BigTorusPolynomial &a, const BigTorusPolynomial &b);

/** @brief substraction */
void sub(BigTorusPolynomial &reps, const BigTorusPolynomial &a, const BigTorusPolynomial &b);


/** @brief out=-in*/
void neg(BigTorusPolynomial &out, const BigTorusPolynomial &in);

/** @brief out=in*/
void copy(BigTorusPolynomial &out, const BigTorusPolynomial &in);

/** @brief int-polynomial torus-polynomial external product  */
void naive_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, UINT64 out_limb_prec);


/** @brief int-polynomial torus-polynomial external product using fft  */
void fft_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, const UINT64 bits_a,
                          UINT64 out_limb_prec);

/** @brief int-polynomial torus-polynomial external product using fft  */
void fft_external_product(BigComplex *out, int64_t *a, BigComplex *in, const int64_t N, const UINT64 bits_a);

/** @brief polynomial-polynomial torus-polynomial internal product using fft  */
void fft_internal_product(BigTorusPolynomial &out, const BigTorusPolynomial &a, const BigTorusPolynomial &b,
                          UINT64 out_limb_prec);

/** @brief int-polynomial torus-polynomial external product using fft  (TODO: a voir comment definir a)*/
void fft_semi_external_product(BigTorusPolynomial &out, const BigComplex *ca, const BigTorusPolynomial &b,
                               const UINT64 bits_a, UINT64 out_limb_prec);

/** @brief polynomial-polynomial torus-polynomial internal product using fft  */
void iFFT(BigComplex *out, const BigTorusPolynomial &a, const UINT64 fft_limbs);

/** @brief polynomial-polynomial torus-polynomial internal product using fft  */
void iFFT(BigComplex *out, const BigTorusPolynomial &a);

/** copy msb bits of each coefficient of b*2^lshift in an array of BigReal */
void precise_conv_toBigReal(BigReal *array, const BigTorusPolynomial &b, int64_t lshift, int64_t msbToKeep);


/** compute  out= array* 2^left_shift mod 1  */
void BigRealPoly_shift_toBigTorus(BigTorusPolynomial &out, BigReal *array, int64_t left_shift);

#endif //FHE_BIGTORUSPOLYNOMIAL_H
