#ifndef FHE_TRLWE_H
#define FHE_TRLWE_H


#include "BigTorusPolynomial.h"
#include "TLwe.h"

class TRLweParams: public TLweParams {
public:
    TRLweParams(const UINT64 N, const BigTorusParams &fixp_params);
    ~TRLweParams() = default;
};

class TRLwe {
public:
    const TRLweParams& params;
    BigTorusPolynomial* const a;

    TRLwe(const TRLweParams& params);

    ~TRLwe();
};


struct pubKsKey {
    TRLwe** kskey;
    BigTorus bitDecomp_in_offset; // sum Bg/2 Bg^i
    __int128 bitDecomp_out_offset; // -Bg/2
    UINT64 BgBits;
    UINT64 l_dec;
    UINT64 out_prec_limbs;
};

void pubKS(TRLwe& out, TLwe& in, pubKsKey& ks);

/**
 * @brief create a trivial ciphertext of a constant
 */
void trivial(TRLwe& out, const BigTorusRef& in, const UINT64 out_limb_prec);

/** @brief get the j-th coef of the base 2^128 decomposition between [0 and 2^128[
 * the output can be signed, we don't care */
__int128 bitdecomp_coef128(const BigTorusRef& tmpDec, UINT64 j, const UINT64 limb_prec);


/** @brief out = out - aij * in
 *  WARNING: input limb precision = out_limb_prec + 2 */
void subMul(TRLwe& out, __int128 aij, const TRLwe& in, const UINT64 out_limb_prec);


#endif //FHE_TRLWE_H
