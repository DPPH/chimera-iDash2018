#ifndef FHE_TRGSW_H
#define FHE_TRGSW_H


#include "BigTorusPolynomial.h"
#include "TRLwe.h"

class TRGSWParams : public TRLweParams {
public:
    static const UINT64 Bgbits = 32;
    static const UINT64 ell = 3;

    TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params);

    ~TRGSWParams() = default;
};

class TRGSW {
public:
    const TRGSWParams &params;
    BigComplex *(a[2][3][2]);

    TRGSW(const TRGSWParams &params);

    ~TRGSW();
};


void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits);

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 bits_a, UINT64 alpha_bits);

#endif //FHE_TRGSW_H
