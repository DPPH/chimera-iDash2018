#ifndef FHE_TRGSW_H
#define FHE_TRGSW_H


#include "BigTorusPolynomial.h"
#include "TRLwe.h"

class TRGSWParams : public TRLweParams {
public:
    static const UINT64 Bgbits = 32;
    static const UINT64 max_ell = 3;

    TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params);

    ~TRGSWParams() = default;
};

class TRGSW {
public:
    const TRGSWParams &params;
    BigComplex *(a[2][3][2]);
    UINT64 bits_a;     //plaintext norm 1 bits
    UINT64 fft_nlimbs; //fft limbs in all BigComplex
    UINT64 ell;        //actual decomposition length

    TRGSW(const TRGSWParams &params);

    ~TRGSW();
};

void intPoly_encrypt(TRGSW &reps, const int64_t *plaintext, const TLweKey &key, UINT64 alpha_bits);

void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits);

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 out_alpha_bits);

#endif //FHE_TRGSW_H
